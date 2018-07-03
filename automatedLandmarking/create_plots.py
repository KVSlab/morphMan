from common import *
from os import path, listdir
from matplotlib.pyplot import *

def prity_title(text):
    text = text.replace("_", " ")
    text = text.capitalize()
    text = text.replace("min", "minimum")
    text = text.replace("max", "maximum")
    return text

def plot_area(data):
    key = ["area", "derivative", "length_min_max", "local_max", "circleness",
            "number_of_max"]
    data_area = {}
    map = {"R": 2, "None": 0, "U": 1}
    color = {0: "bo", 1: "g^", 2:"ys"}

    folders = data.keys()
    for k in data[folders[0]].keys():
        for i in range(len(key)):
            if key[i] in k:
                data_area[k] = {"R": [], "U": [], "None": []}

    for folder in folders:
        for k in data_area.keys():
            if "C0" in folder:
                rupture = data[folder]["ruptureStatus"]
                rupture = rupture if "function" not in rupture else "R"
            else:
                rupture = "None"

            if k != "local_max_stepest_disent":
                data_area[k][rupture].append(data[folder][k])
            else:
                data_area[k][rupture].append(data[folder][k][0])

    for k, v in data_area.iteritems():
        figure()
        # Set up plot
        xticks([0, 1, 2], ("Normal", "Unruptured", "Ruptured"))
        title(prity_title(k))

        results = []
        for result in data_area[k].values():
            results += result

        # Set range of the plot
        xlim([-0.5, 2.5])
        max_ = (np.asarray(results)*1.05).max()
        min_ = (np.asarray(results)*0.40).min()
        #ylim([min_, max_])

        for type, value in data_area[k].iteritems():
            # Get case type
            status = map[type]

            # Plot
            for v in value:
                plot([status], v, color[status])
                hold("on")

        # Plot mean and standard deviation
        for key in ["None", "U", "R"]:
            tmp = np.asarray(data_area[k][key])
            tmp_mean = np.mean(tmp)
            SD = np.std(tmp)
            status = map[key]
            fmt = color[status][1]
            color_ = color[status][0]

            errorbar([status], [tmp_mean], color="k", fmt=fmt, linewidth=2, 
                     yerr=[[SD], [SD]])
            hold("on")

        # Show / save
        savefig(k + ".png")
        #show()


if __name__ == "__main__":
    #basedir = "."
    basedir = "/home/aslak/master/src/aneurysms"
    folders = listdir(basedir)
    data = {}
    for folder in folders:
        if "P0" in folder or "C0" in folder:
            data[folder] = getParameters(path.join(basedir, folder))

    basedir = path.join(basedir, "new_cases")
    folders = listdir(basedir)
    for folder in folders:
        if "P0" in folder or "C0" in folder:
            if folder in ["C0023", "C0099", "C0057b", "C0093"]: continue
            data[folder] = getParameters(path.join(basedir, folder))

    plot_area(data)
