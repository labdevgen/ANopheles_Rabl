from Expected_calculator import get_dumpPath, get_contacts_using_juicer_dump
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os,pickle

hic_files = ["https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AalbS2/AalbS2_V4_30.hic",
            "https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AatrE3/AatrE3_V4_30.hic",
            "https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AsteI2/AsteI2_V4_30.hic",
             "https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AcolNg/AcolNg_V4.hic",
             "https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AmerR4/AmerR4_V4_30.hic",
            "https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AmerR4/AmerR4A_V4_30.hic"]

resolution = 50000
juicer_path = os.path.abspath("../Anopheles_Ps/Juicer/juicer_tools_1.19.02.jar")
wings = [("2R","2L"),("3R","3L")]

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
plt.gcf().set_size_inches(14,10)

fout = open("results.txt","w")
for hic_file in hic_files:
    curr_rabl = []
    curr_others = []
    for ind,(chr1,chr2) in enumerate(wings):
        contacts = dump_contacts(hic_file, resolution, chr1, chr2, juicer_path)
        if not (chr1.endswith("R") and chr2.endswith("L")):
            # Note: if you change order of chrms
            # you would need to change coordinate system of contacts differentely
            raise

        max_L_size = max(contacts.st1.values)

        contacts.dropna(inplace=True)

        contacts.st1 = contacts.apply(lambda x: max_L_size - x)
        contacts["count"] = contacts["count"].apply(np.log)
        contacts["rabl_dist"] = (contacts.st1 - contacts.st2).abs()
        contacts["cent_dist"] = (contacts.st1 + contacts.st2)

        rabl = contacts.query("rabl_dist <= 3*@resolution")["count"].values
        non_rabl = contacts.query("rabl_dist >= 6*@resolution")["count"].values
        curr_rabl += rabl.tolist()
        curr_others += non_rabl.tolist()

        av_rabl = np.exp(np.average(rabl))
        av_nonrabl = np.exp(np.average(non_rabl))
        fout.write(hic_file + "\t" + chr1+"-"+chr2+"\t"+str(av_rabl/av_nonrabl)+"\n")

    data += [curr_rabl,curr_others]
    labels += ["Rabl","Other"]

plt.axhline(y=0, ls=":")
plt.vlines(x=(np.arange(1, len(hic_files)))*2+0.5, ymin=-3, ymax=3)
plt.boxplot(data,labels=labels,showfliers=False)
plt.xticks(rotation=45, fontsize = 12)
plt.ylabel("Log(normed contacts count)", fontsize = 14)
for ind,val in enumerate(["An. albimanus","An. atroparvus",
                          "An. stephensi","An. coluzzii",
                          "An. merus\n(embryo)","An. merus\n(adult)"]):
    plt.annotate(val,(ind*2+1.5,-2.5), style="italic", size=14,
                 horizontalalignment="center",
                 verticalalignment="top")
plt.savefig("Fig_rabl.png",dpi=150)
