
import argparse as ap
import matplotlib.pyplot as plt

def main():
    # plot the minimum of n categories

    parser = ap.ArgumentParser(description="options")

    parser.add_argument("--minima","-m", type=float, default=[], nargs="*" )
    parser.add_argument("--output","-o", type=str, default="")

    args = parser.parse_args()

    x_vals = [i for i in range(2, 2+len(args.minima))]

    fig, ax = plt.subplots(1,1)

    ax.plot(x_vals, args.minima, label="Category Minima", marker='.')

    ax.set_xlabel("Number of Categories")
    plt.xticks([ 2, 3, 4, 5, 6])

    ax.set_ylabel("Minimum [A.U.]")
    ax.legend(loc='best')

    fig.savefig(f'optNumCats_{args.output}.png', bbox_inches="tight")
    

if __name__ == "__main__":
    main()

