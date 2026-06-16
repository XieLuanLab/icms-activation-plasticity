import matplotlib as mpl
import seaborn as sns
import matplotlib.pyplot as plt

PALETTE = sns.color_palette("deep")


def apply_global_style():
    mpl.rcParams.update({
        # --- Fonts ---
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial'],   # requires Arial installed

        # --- Sizes (final print target ~5–7 pt ticks) ---
        'font.size': 7,
        'axes.titlesize': 7,
        'figure.titlesize': 7,
        'legend.fontsize': 7,
        'axes.labelsize': 7,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,

        # --- Lines & ticks (thinner) ---
        'axes.linewidth': 0.32,      # ~0.3 pt spines
        'lines.linewidth': 0.8,      # data lines
        'lines.markersize': 2.5,
        'lines.markeredgewidth': 0,
        'errorbar.capsize': 5,

        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.size': 2.5,     # length in points
        'ytick.major.size': 2.5,
        'xtick.major.width': 0.32,
        'ytick.major.width': 0.32,

        # --- Spines & legend ---
        'axes.spines.top': False,
        'axes.spines.right': False,
        'legend.frameon': False,

        # --- Export: keep text editable ---
        'svg.fonttype': 'none',

        # --- Save defaults ---
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.01,
    })


if __name__ == "__main__":
    apply_global_style()

    fig, ax = plt.subplots(figsize=(2.5, 2.5))  # ~63.5 mm square
    ax.plot([0, 1, 2], [0, 1, 4], label="demo")
    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.set_title("Hello world!")
    ax.legend()

    # Optional belt-and-suspenders for patches:
    # ax.set_facecolor('none')
    # fig.patch.set_alpha(0)
    # ax.patch.set_linewidth(0)
    # ax.patch.set_edgecolor('none')

    # Save (rcParams already set bbox/pad)
    plt.savefig("C:/Users/robin/Downloads/figure.svg", format="svg")
    # (Optional backup)
    # plt.savefig("C:/Users/robin/Downloads/figure.pdf", format="pdf")
