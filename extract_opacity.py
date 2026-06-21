import h5py
import numpy as np

def extract_to_text(h5_path, out_path, z_val="0.020000"):
    with h5py.File(h5_path, 'r') as f:
        logTs = f['logTs'][:]
        logRs = f['logRs'][:]
        Xs = f[z_val]['Xs'][:]
        # Shape: (logT, logR, X, fCO, fC, fN)
        # We take the base mixture (index 0 for the last 3 axes)
        kap = f[z_val]['kap'][:, :, :, 0, 0, 0]

        with open(out_path, 'w') as out:
            out.write(f"AESOPUS / AGSS09 Opacity Table (Z = {float(z_val)})\n")
            out.write("====================================================\n\n")
            
            for i, x in enumerate(Xs):
                out.write(f"Table for X = {x:.3f}\n")
                out.write(f"logT / logR | " + " ".join([f"{r:6.3f}" for r in logRs]) + "\n")
                out.write("-" * 150 + "\n")
                
                for j, t in enumerate(logTs):
                    row = [f"{t:5.3f}     |"]
                    for k in range(len(logRs)):
                        val = kap[j, k, i]
                        row.append(f"{val:6.3f}")
                    out.write(" ".join(row) + "\n")
                out.write("\n")

extract_to_text('/Users/alec/Downloads/AESOPUS_AGSS09.h5', 'opacity_AGSS09_z0.02.txt')
print("Successfully extracted opacity table to opacity_AGSS09_z0.02.txt")
