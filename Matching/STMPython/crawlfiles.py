import os
for root, dirs, files in os.walk("/Xnfs/site/lagrangian/PTV/"):
    for file in files:
        if file.endswith(“.dat”):
             print(os.path.join(root, file))
