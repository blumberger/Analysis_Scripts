import os

if not os.path.isdir("snapshots/old"):
    os.makedirs("snapshots/old")

for i in range(int(count)): 
    i = str(i).zfill(5) 
    f1 = f"snapshots/temp_frame_{i}.png" 
    f2 = f"snapshots/ang_dist_frame_{i}.png" 
    os.system(f"convert +append {f1} {f2} snapshots/frame_{i}.jpg")
    os.rename(f1, f"snapshots/old/temp_frame_{i}.png")
    os.rename(f2, f"snapshots/old/ang_dist_frame_{i}.png")

