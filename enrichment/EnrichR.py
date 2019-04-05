import argparse

from enrichment.regulon_enrichment import *


def main():

    parser = argparse.ArgumentParser(
        "Generate regulon enrichment scores for an expression dataset using Pathway Commons features ")
    parser.add_argument('--expr_f', metavar = '', help='tab-separated file of expression data of\
     shape = [n_features, n_samples]')
    parser.add_argument('--cohort', metavar = '', help='name of cohort to associate with compiled regulon and\
     enrichment scorestab-separated file of expression data of shape = [n_features, n_samples]')
    parser.add_argument('--norm_type', metavar = '', help='Scaling method to implement:\
     robust | standard | minmax | quant', default='robust')
    parser.add_argument("--feature", action = "store_true", help = "Scale across features", default = True)
    parser.add_argument("--sample", action = "store_true", help = "Scale across features and samples", default = False)
    parser.add_argument("--thresh_filter", metavar = '', help = "Filter non-expressed features in dataset\
     before scaling", default = 0.4)
    parser.add_argument("--scale", action = "store_true", help = "optional arg to avoid scaling dataset if data set has\
     been normalized prior to analysis", default = True)
    parser.add_argument('--regulon_size', type=int, metavar = '', help='Minimum downstream regulon members',
                        default = 15)

    args = parser.parse_args()
    print(args)
    generate_enrichment_scores(expr_f = args.expr_f, cohort = args.cohort, norm_type = args.norm_type,
                               feature = args.feature, sample = args.sample, thresh_filter = args.thresh_filter,
                               scale = args.scale, regulon_size = args.regulon_size)


if __name__ == '__main__':
    main()
