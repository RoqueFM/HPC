#! /usr/bin/env python3

import sys
import numpy as np



import matplotlib.pyplot as plt
import lib.plot_config as pc

if len(sys.argv) < 2:
    print("")
    print("Usage:")
    print("")
    print(f" {sys.argv[0]} [output_file.csv]")
    print("")
    sys.exit(1)

fig, ax = pc.setup(scale=1.5)

ax.set_title(", ".join(sys.argv[1:]))
ax.set_xlabel("N")
ax.set_ylabel("GFlops/s")
ax.set_xscale("log")

ps = pc.PlotStyles()


for csvfile in sys.argv[1:]:
    g = csvfile
    g = g.replace("output_run_", "")
    g = g.replace("mmul_", "")
    g = g.replace(".csv", "")

    rawdata = np.genfromtxt(csvfile, delimiter='\t', dtype=str)

    data = np.array(rawdata[1:,1:], dtype=float)
    labels = rawdata[1:,0]
    ticks = np.array(rawdata[0,1:], dtype=float)

    for i in range(len(labels)):
        plot_style = ps.getNextStyle()
        label = labels[i]
        label = label.replace("mm_mul_", "")
        label = label.replace("matrix_matrix_mul_", "")

        if len(sys.argv) >= 2:
            label = g+" "+label

        plt.plot(ticks, data[i,:], label=label, **plot_style)


ax.legend()

output_file = sys.argv[1]

output_file = output_file.replace('.csv', '.pdf')

if output_file != sys.argv[1]:
    print("Writing output to file: "+output_file)
    pc.savefig(output_file)
    print("Done")

plt.show()

