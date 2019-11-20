from __future__ import with_statement
import os
import numpy as np

if __name__ == "__main__":

    source_dir = os.path.join("/Users","danielsf","physics","newSedLibrary")
    source_dir = os.path.join(source_dir, "starSED", "mlt")
    file_name_list = os.listdir(source_dir)
    target_dir = os.path.join("/Users", "danielsf", "physics", "newSedLibrary")
    target_dir = os.path.join(target_dir, "starSED", "phoSimMLT")


    dtype = np.dtype([("wavelen", np.float), ("flux", np.float)])

    for file_name in file_name_list:
        data = np.genfromtxt(os.path.join(source_dir, file_name), dtype=dtype)
        target_name = os.path.join(target_dir, file_name).replace('.gz','')
        with open(target_name, "w") as output_file:
            for ww, ff in zip(data["wavelen"][:24000], data["flux"][:24000]):
                output_file.write("%.7e %.7e\n" % (ww, ff))
        os.system("gzip %s" % target_name)
