import os
import gc
import time
import datetime
import psutil
import pandas as pd
import numpy as np
import h5py
from itertools import combinations
from tqdm import tqdm
from multiprocessing import Pool, current_process

# ----------------------------
# Global variables for workers
# ---------------------------
DF_FILE = None       # Each worker will read from this CSV
DF_INDEX_COL = 0     # The first column contains gene names (index)
THRESHOLD = 0        # Threshold for expression data
CHECKPOINT_DIR = "checkpoints"

def init_worker(df_path):
    """
    Initialize worker with file path only - don't load data yet
    """
    global DF_FILE
    DF_FILE = df_path

def check_memory_usage():
    """
    Check memory usage and pause if it's too high
    """
    if psutil.virtual_memory().percent > 80:  # If over 80% memory used
        process_name = current_process().name
        print(f"[Worker {process_name}] High memory usage detected ({psutil.virtual_memory().percent}%). Pausing for 60 seconds.")
        time.sleep(60)  # Wait for GC to run
        gc.collect()    # Force collection
        print(f"[Worker {process_name}] Resuming after memory cleanup. Current usage: {psutil.virtual_memory().percent}%")

def process_chunk(args):
    """
    Process a chunk of genes with memory optimization
    """
    genes_chunk, chunk_idx = args
    checkpoint_file = f"{CHECKPOINT_DIR}/chunk_{chunk_idx}.npy"
   
    if os.path.exists(checkpoint_file):
        print(f"[Worker {current_process().name}] Found existing checkpoint for chunk {chunk_idx}, skipping.")
        return None

    print(f"[Worker {current_process().name}] Processing chunk {chunk_idx} with {len(genes_chunk)} genes")
    check_memory_usage()  # Check memory before starting a new chunk

    # Create a local square matrix for this chunk
    chunk_size = len(genes_chunk)
    local_mat = np.zeros((chunk_size, chunk_size), dtype=np.int32)
   
    # Map each gene in the chunk to its index in the matrix
    gene_to_idx = {gene: i for i, gene in enumerate(genes_chunk)}
    genes_set = set(genes_chunk)
    
    # First read the header to get column names
    df_header = pd.read_csv(DF_FILE, nrows=0)
    all_columns = df_header.columns.tolist()
    
    # Read the CSV in chunks to reduce memory usage
    chunk_size_read = 1000  # Read 1000 samples at a time
    total_cols_processed = 0
    
    # Use all columns instead of a lambda function
    for df_chunk in pd.read_csv(DF_FILE, chunksize=chunk_size_read, index_col=0, dtype={0: str}):
  
        # Check memory before processing each dataframe chunk
        check_memory_usage()
        
        # Process only the genes in our chunk
        df_filtered = df_chunk.loc[df_chunk.index.intersection(genes_set)]
        
        if df_filtered.empty:
            continue
        
        # Iterate over all samples (columns) in this chunk of the DataFrame
        for col in df_filtered.columns:
            # Find all genes not expressed in this sample (<= threshold)
            not_expressed = df_filtered.index[df_filtered[col] <= THRESHOLD].tolist()
            
            # Update diagonal counts
            for g in not_expressed:
                idx = gene_to_idx[g]
                local_mat[idx, idx] += 1

            # Update pairwise counts
            for g1, g2 in combinations(not_expressed, 2):
                idx1 = gene_to_idx[g1]
                idx2 = gene_to_idx[g2]
                local_mat[idx1, idx2] += 1
                local_mat[idx2, idx1] += 1
        
        total_cols_processed += len(df_chunk.columns)
        if total_cols_processed % 1000 == 0:
            print(f"[Worker {current_process().name}] Chunk {chunk_idx}: processed {total_cols_processed} columns")
        
        # Force garbage collection after each dataframe chunk
        del df_filtered
        del df_chunk
        gc.collect()

    # Save the result of this chunk to a checkpoint file
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    np.save(checkpoint_file, local_mat)
    print(f"[Worker {current_process().name}] Saved checkpoint for chunk {chunk_idx}")
    
    del local_mat
    del gene_to_idx
    gc.collect()  # Explicit garbage collection
    return None

def main():
    process_info = psutil.Process()
    print(f"Initial memory usage: {process_info.memory_info().rss / (1024**2):.2f} MB")
    start_time = datetime.datetime.now()
    print(f"Starting processing at: {start_time}")

    # -----------------
    # Configuration
    # -----------------
    expression_csv = "/home/user-kp/anugreha/human/test/erythrocyte_test.csv"
    gene_ids_file = "/home/user-kp/anugreha/human/test/symbols_test.csv"
    hdf5_file = "/home/user-kp/anugreha/human/test/final_matrix_erythrocyte_test.h5"
        
    chunk_size = 100      # Smaller chunk size to reduce memory usage
    num_processes = 4     # Adjust based on your system's capabilities
   
    print("Loading gene IDs...")
    gene_ids = pd.read_csv(gene_ids_file)['Gene'].tolist()
    total_genes = len(gene_ids)
    print(f"Loaded {total_genes} genes from {gene_ids_file}")

    print("Checking expression data dimensions...")
    # Read the first row to get column count
    df_header = pd.read_csv(expression_csv, nrows=0)
    total_samples = len(df_header.columns) - 1  # -1 for the index column
    
    print(f"Expression data dimensions (approx): {total_genes} genes x {total_samples} samples")
    print(f"Current memory usage: {psutil.virtual_memory().percent}% of total RAM")
    
    # Check memory before creating chunks
    check_memory_usage()

    num_chunks = (total_genes + chunk_size - 1) // chunk_size
    last_chunk_size = total_genes - (num_chunks - 1) * chunk_size

    print(f"Processing will happen in {num_chunks} chunks.")
    print(f"Regular chunk size: {chunk_size}")
    print(f"Last chunk size:   {last_chunk_size}\n")

    # Prepare chunk tasks
    tasks = []
    for start_idx in range(0, total_genes, chunk_size):
        chunk_idx = start_idx // chunk_size
        end_idx = min(start_idx + chunk_size, total_genes)
        genes_chunk = gene_ids[start_idx:end_idx]
        print(f"Chunk {chunk_idx} size: {len(genes_chunk)} genes")
        tasks.append((genes_chunk, chunk_idx))

    # Create checkpoint directory if it doesn't exist
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)

    # ---------------
    # Multiprocessing
    # ---------------
    print("\n[Main] Starting multiprocessing Pool ...")
    with Pool(processes=num_processes,
              initializer=init_worker,
              initargs=(expression_csv,)) as pool:
        # Process chunks in order of increasing size for better resource management
        for _ in tqdm(pool.imap_unordered(process_chunk, tasks),
                      total=len(tasks), desc="Processing chunks"):
            # Check memory after each chunk completes
            if psutil.virtual_memory().percent > 80:
                print(f"[Main] High memory usage detected ({psutil.virtual_memory().percent}%). Pausing for 60 seconds.")
                time.sleep(60)
                gc.collect()
                print(f"[Main] Resuming after memory cleanup. Current usage: {psutil.virtual_memory().percent}%")
   
    # -------------------------------------------
    # Combine all checkpoint files into one HDF5 incrementally
    # -------------------------------------------
    print("\nCombining chunk results into final HDF5 matrix incrementally...")
    with h5py.File(hdf5_file, "w") as f:
        # Create the dataset with chunked storage for efficient partial I/O
        dset = f.create_dataset("final_matrix",
                                (total_genes, total_genes),
                                dtype="int32",
                                chunks=(min(1000, total_genes), min(1000, total_genes)),
                                compression="gzip",
                                compression_opts=4)
       
        for chunk_idx, (genes_chunk, _) in enumerate(tasks):
            # Check memory before loading each checkpoint
            check_memory_usage()
            
            checkpoint_file = f"{CHECKPOINT_DIR}/chunk_{chunk_idx}.npy"
            if not os.path.exists(checkpoint_file):
                print(f"Warning: Missing checkpoint file for chunk {chunk_idx}")
                continue
                
            # Calculate the indices in the global matrix
            start_idx = chunk_idx * chunk_size
            end_idx = min(start_idx + len(genes_chunk), total_genes)
            
            # Load the checkpoint file
            chunk_mat = np.load(checkpoint_file)
            
            # Write to the appropriate location in the HDF5 file
            dset[start_idx:end_idx, start_idx:end_idx] = chunk_mat
            
            # Clean up
            del chunk_mat
            gc.collect()
            
            # Report progress
            print(f"Added chunk {chunk_idx} to HDF5 file ({start_idx}:{end_idx})")

    print("Final matrix saved in HDF5 format at:", hdf5_file)

    end_time = datetime.datetime.now()
    runtime = end_time - start_time
    print("\nProcessing completed.")
    print(f"Start time:   {start_time}")
    print(f"End time:     {end_time}")
    print(f"Total runtime: {runtime}")

if __name__ == "__main__":
    main()
