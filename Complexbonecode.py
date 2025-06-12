import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import deque
import time


# Step 1: Read and process strain tensor data
data = []
with open('3D_strain_output_v1.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('%'):
            continue
        tokens = line.split()
        if len(tokens) >= 9:
            try:
                x, y, z = map(float, tokens[:3])
                eXX, eXY, eXZ, eYY, eYZ, eZZ = map(float, tokens[3:9])
            except (ValueError, TypeError):
                continue
            
            if np.isnan([eXX, eXY, eXZ, eYY, eYZ, eZZ]).any():
                continue
            
            data.append([x, y, z, eXX, eXY, eXZ, eYY, eYZ, eZZ])

start_time = time.time()#initial time elapsed



data = np.array(data)
X, Y, Z = data[:, 0], data[:, 1], data[:, 2]
eXX, eXY, eXZ, eYY, eYZ, eZZ = data[:, 3], data[:, 4], data[:, 5], data[:, 6], data[:, 7], data[:, 8]

# Step 2: Create structured grid
x_vals, y_vals, z_vals = np.unique(X), np.unique(Y), np.unique(Z)
nx, ny, nz = len(x_vals), len(y_vals), len(z_vals)

x_to_i = {x: i for i, x in enumerate(x_vals)}
y_to_j = {y: j for j, y in enumerate(y_vals)}
z_to_k = {z: k for k, z in enumerate(z_vals)}

# Step 3: Store strain tensors and mark valid points
A = np.zeros((nx, ny, nz, 3, 3))
valid = np.zeros((nx, ny, nz), dtype=bool)
for row in data:
    i = x_to_i[row[0]]
    j = y_to_j[row[1]]
    k = z_to_k[row[2]]
    A[i, j, k] = [[row[3], row[4], row[5]],
                  [row[4], row[6], row[7]],
                  [row[5], row[7], row[8]]]
    valid[i, j, k] = True

# Step 4: Compute eigenvectors for valid points
Q = np.zeros((nx, ny, nz, 3, 3))
skipped_count = 0

for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            if valid[i, j, k]:
                try:
                    vals, vecs = np.linalg.eigh(A[i, j, k])
                    idx = np.argsort(vals)[::-1]
                    Q[i, j, k] = vecs[:, idx]
                except np.linalg.LinAlgError:
                    valid[i, j, k] = False
                    skipped_count += 1

# Step 5: Propagate eigenvector consistency using grid order
qupdated = np.zeros((nx, ny, nz), dtype=bool)
neighbors = [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)]

def best_match(ref_vec, options):
    dots = [np.dot(ref_vec, opt) for opt in options]
    return options[np.argmax(dots)]

#Queue for the BFS algorithm
queue = deque()
if valid[0, 0, 0]:
    qupdated[0, 0, 0] = True
    queue.append((0, 0, 0))

while queue:
    i, j, k = queue.popleft()
    q_ref = Q[i, j, k]
    for di, dj, dk in neighbors:
        ni, nj, nk = i + di, j + dj, k + dk
        if 0 <= ni < nx and 0 <= nj < ny and 0 <= nk < nz:
            if valid[ni, nj, nk] and not qupdated[ni, nj, nk]:
                current_q = Q[ni, nj, nk]
                q1_options = [
                    current_q[:, 0], -current_q[:, 0],
                    current_q[:, 1], -current_q[:, 1],
                    current_q[:, 2], -current_q[:, 2]
                ]
                q2_options = [
                    current_q[:, 1], -current_q[:, 1],
                    current_q[:, 0], -current_q[:, 0],
                    current_q[:, 2], -current_q[:, 2]
                ]
                #Calculate the most consistent vector for all q components
                best_q1 = best_match(q_ref[:, 0], q1_options)
                best_q2 = best_match(q_ref[:, 1], q2_options)
                best_q3 = np.cross(best_q1, best_q2)
                
                Q[ni, nj, nk, :, 0] = best_q1
                Q[ni, nj, nk, :, 1] = best_q2
                Q[ni, nj, nk, :, 2] = best_q3
                qupdated[ni, nj, nk] = True
                queue.append((ni, nj, nk))

end_time = time.time()
execution_time = end_time - start_time
print(execution_time)
# Flip q3 of the very first point (after all propagation)
if valid[0, 0, 0]:
    Q[0, 0, 0, :, 2] *= -1  # Flip direction
    # Optional: Recompute q3 to ensure orthogonality
    Q[0, 0, 0, :, 2] = np.cross(Q[0, 0, 0, :, 0], Q[0, 0, 0, :, 1])

for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            if valid[i, j, k]:
                Q[i, j, k, :, 0] *= -1  # Flip q1
                Q[i, j, k, :, 1] *= -1  # Flip q2


# Save result
output_path = "eigenvectors_bfs_output1.1.txt"
with open(output_path, "w") as f:
    f.write("X Y Z Q00 Q01 Q02 Q10 Q11 Q12 Q20 Q21 Q22\n")
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if valid[i, j, k]:
                    x, y, z = x_vals[i], y_vals[j], z_vals[k]
                    q_values = Q[i, j, k].flatten()
                    f.write(f"{x} {y} {z} " + " ".join(map(str, q_values)) + "\n")

print(f"Eigenvectors saved to {output_path}")

output_path, skipped_count
# Step 7: Visualization with valid points
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
sample_rate = max(1, nx//20)

for i in range(0, nx, sample_rate):
    for j in range(0, ny, sample_rate):
        for k in range(0, nz, sample_rate):
            if not valid[i, j, k]:
                continue
            x, y, z = x_vals[i], y_vals[j], z_vals[k]
            q1, q2, q3 = Q[i, j, k, :, 0], Q[i, j, k, :, 1], Q[i, j, k, :, 2]
            ax.quiver(x, y, z, q1[0], q1[1], q1[2], color='b', length=0.02)
            ax.quiver(x, y, z, q2[0], q2[1], q2[2], color='g', length=0.02)
            ax.quiver(x, y, z, q3[0], q3[1], q3[2], color='r', length=0.02)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title('Consistent Eigenvectors')
plt.show()