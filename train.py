import numpy as np
import pandas as pd

# Define which columns to use for features and target based on the project description
# machine.data columns (0-indexed):
# 2: MYCT
# 3: MMIN
# 4: MMAX
# 5: CACH
# 6: CHMIN (skipped in the project model)
# 7: CHMAX
# 8: PRP (target)
feature_indices = [2, 3, 4, 5, 7] # MYCT, MMIN, MMAX, CACH, CHMAX
target_index = 8                   # PRP

raw_data_records = []
with open('data/machine.data', 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split(',')
        
        # Ensure the line has enough columns before trying to access them
        if len(parts) < max(feature_indices + [target_index]) + 1:
            print(f"Skipping line due to insufficient columns: {line}")
            continue
            
        try:
            record_features = [float(parts[i]) for i in feature_indices]
            record_target = float(parts[target_index])
            raw_data_records.append(record_features + [record_target])
        except ValueError:
            print(f"Skipping line due to invalid numeric data: {line}")
            continue

if not raw_data_records:
    print("Error: No data records were successfully read. Exiting.")
    exit(1)

data_np = np.array(raw_data_records)
n_rows, n_cols_data = data_np.shape  # n_cols_data will be num_features + 1
num_features = len(feature_indices)

print(f"Successfully read {n_rows} records.")
print(f"Number of features: {num_features}")
print(f"Data shape (rows, features + target): {n_rows}, {n_cols_data}")

# Sequential train/test split (80% train, 20% test)
train_size = int(n_rows * 0.8)
test_size = n_rows - train_size

print(f"Total samples: {n_rows}")
print(f"Training samples: {train_size}")
print(f"Testing samples: {test_size}")

if train_size == 0 or test_size == 0:
    print("Error: Train or test set size is 0. Adjust data or split ratio.")
    exit(1)

# No shuffling, take data sequentially
X_train = data_np[:train_size, :num_features]
y_train = data_np[:train_size, num_features]
X_test = data_np[train_size:, :num_features]
y_test = data_np[train_size:, num_features]

print(f"X_train shape: {X_train.shape}")
print(f"y_train shape: {y_train.shape}")
print(f"X_test shape: {X_test.shape}")
print(f"y_test shape: {y_test.shape}")

# Build normal equations: A_ls * beta = b_ls
# A_ls = X_train^T * X_train
# b_ls = X_train^T * y_train
A_ls = X_train.T.dot(X_train)
b_ls = X_train.T.dot(y_train)

print("\nMatrix A_ls (X_train^T * X_train):")
print(A_ls)
print("\nVector b_ls (X_train^T * y_train):")
print(b_ls)

# Solve for model parameters (beta)
try:
    parameters_beta = np.linalg.solve(A_ls, b_ls)
except np.linalg.LinAlgError:
    print("Error: Singular matrix encountered. Cannot solve the linear system.")
    print("This might happen if features are perfectly collinear or A_ls is ill-conditioned.")
    print("Consider using np.linalg.lstsq for a least-squares solution if this is expected.")
    exit(1)


print("\nModel Parameters (beta from Python):")
feature_names_ordered = ["MYCT", "MMIN", "MMAX", "CACH", "CHMAX"]
for i in range(num_features):
    print(f"x{i+1} ({feature_names_ordered[i]}): {parameters_beta[i]:.4f}")

# Make predictions on the test set
y_pred_test = X_test.dot(parameters_beta)

# Calculate Root Mean Square Error (RMSE) on the test set
mse_test = np.mean((y_pred_test - y_test)**2)
rmse_test = np.sqrt(mse_test)
print(f"\nRoot Mean Square Error (RMSE) on Testing Set: {rmse_test:.4f}")

# Make predictions on the training set
y_pred_train = X_train.dot(parameters_beta)

# Calculate Root Mean Square Error (RMSE) on the training set
mse_train = np.mean((y_pred_train - y_train)**2)
rmse_train = np.sqrt(mse_train)
print(f"Root Mean Square Error (RMSE) on Training Set: {rmse_train:.4f}")

print("\nPython script finished.")
