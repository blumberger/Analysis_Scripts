import os

if 'count' not in globals():
   count = 10000

snapshots_folder = "snapshots"

if not os.path.isdir(f"{snapshots_folder}/old"):
    os.makedirs(f"{snapshots_folder}/old")

for i in range(int(count)): 
    i = str(i).zfill(5) 

    if not os.path.isfile(f"{snapshots_folder}/temp_frame_{i}.png"):
      break

    f1 = f"{snapshots_folder}/temp_frame_{i}.png" 
    f2 = f"{snapshots_folder}/ang_dist_frame_{i}.png" 
    f3 = f"{snapshots_folder}/RDF_frame_{i}.png" 
    out_file = f"{snapshots_folder}/frame_{i}.jpg"

    os.system(f"convert +append {f2} {f3} {out_file}")
    os.system(f"convert -append {out_file} {f1} {out_file}")

    os.rename(f1, f"{snapshots_folder}/old/temp_frame_{i}.png")
    os.rename(f2, f"{snapshots_folder}/old/ang_dist_frame_{i}.png")
    os.rename(f3, f"{snapshots_folder}/old/RDF_frame_{i}.png")

print("All frames created!")
