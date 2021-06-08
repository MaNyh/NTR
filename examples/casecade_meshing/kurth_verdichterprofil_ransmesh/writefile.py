f="PointCloudBladeGeometry.txt"

lines = []
with open(f,"r") as fobj:
    lines = fobj.readlines()

new = []
for l in lines:
    ln = l.replace("\n"," 0 \n")
    new.append(ln)

with open("new.txt","w") as fobj:
    fobj.writelines(new)
