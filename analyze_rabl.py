from Expected_calculator import get_dumpPath, get_contacts_using_juicer_dump
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os,pickle

hic_files = ["https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AalbS2/AalbS2_V3_30.hic",
            "https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AatrE3/AatrE3_V3_30.hic",
            "https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AsteI2/AsteI2_V3_30.hic",
             "https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AcolNg/AcolNg_V3_30.hic",
             "https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AmerR4/AmerR4_V3_30.hic",
            "https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AmerR4/AmerR4_V3A_30.hic"]

resolution = 50000
juicer_path = "../Anopheles_Ps/Juicer/juicer_tools_1.19.02.jar"
wings = [("2L","2R"),("3L","3R")]

def dump_contacts(hic_file, resolution, chr1, chr2, juicer_path):
    # dump contacts from chromosome to temp file, then read if to Dframe
    contacts_file = get_dumpPath(hic_file, resolution, chr1,chr2)
    if not os.path.isfile(contacts_file):
        contacts = get_contacts_using_juicer_dump(juicer_path,hic_file,chr1,resolution,chr2=chr2,
                                                  datatype="oe")
        pickle.dump(contacts,open(contacts_file,"wb"))
    return pickle.load(open(contacts_file,"rb"))

data = []
labels = []
for hic_file in hic_files:
    for ind,(chr1,chr2) in enumerate(wings):
        contacts = dump_contacts(hic_file, resolution, chr1, chr2, juicer_path)
        if not (chr1.endswith("L") and chr2.endswith("R")):
            # Note: if you change order of chrms
            # you would need to change coordinate system of contacts differentely
            raise

        if "AcolNg_V3_30" in hic_file:
            print (contacts)
        max_L_size = max(contacts.st1.values)

        contacts.dropna(inplace=True)

        contacts.st1 = contacts.apply(lambda x: max_L_size - x)
        contacts["count"] = contacts["count"].apply(np.log)
        contacts["rabl_dist"] = (contacts.st1 - contacts.st2).abs()
        contacts["cent_dist"] = (contacts.st1 + contacts.st2)

        rabl = contacts.query("rabl_dist <= 3*@resolution")["count"].values
        non_rabl = contacts.query("rabl_dist >= 6*@resolution")["count"].values
        data += [rabl,non_rabl]
        labels += ["Rabl "+chr1+"/"+chr2,"Other "+chr1+"/"+chr2]

plt.axhline(y=0, ls=":")
plt.vlines(x=(np.arange(1, len(hic_files)))*4+0.5, ymin=-3, ymax=3)
plt.boxplot(data,labels=labels,showfliers=False)
plt.xticks(rotation=45)
plt.ylabel("Log(observed/expected)")
plt.show()
