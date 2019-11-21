from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
import matplotlib.pyplot as plt
from decimal import Decimal
from matplotlib.patches import Patch
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from sklearn import linear_model
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold, cross_val_predict, RepeatedKFold, cross_val_score
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
### bbox_to_anchor=(1.,1.),bbox_transform=plt.gcf().transFigure, <- adds legend to top right corner of plot


def plot_enrichment(meta, expr, total_enrichment, regulon, regulator):
    signature = expr.mean() / expr.std()

    control_samps = meta[(meta.condition.isin(['Individual_MNC']))]['RNA_SeqID']
    test_samps = meta[(meta.condition.isin(['ctd2']))]['RNA_SeqID']

    control = total_enrichment.loc[control_samps]
    test = total_enrichment.loc[test_samps]

    activity = test.mean() - control.mean()

    control_expr = expr.loc[control_samps]
    test_expr = expr.loc[test_samps]

    ttest_results_expr = st.ttest_ind(control_expr, test_expr)
    stats_frame_expr = pd.DataFrame(np.array([ttest_results_expr.pvalue, ttest_results_expr.statistic]),
                                    index = ['pvalue', 't-statistic'], columns = control_expr.columns).T

    weighted_ranks = (st.norm.ppf(stats_frame_expr['pvalue'] / 2.0, loc = 0, scale = 1)) * np.sign(
        stats_frame_expr['t-statistic'])
    weighted_ranks.sort_values(inplace = True)
    signature = signature.loc[weighted_ranks.index]

    expression = test_expr.mean() - control_expr.mean()

    ttest_results = st.ttest_ind(control, test)
    stats_frame = pd.DataFrame(np.array([ttest_results.pvalue, ttest_results.statistic]),
                               index = ['pvalue', 't-statistic'], columns = control.columns).T
    stats_frame.sort_values('pvalue', inplace = True)
    heat_frame = pd.DataFrame(np.array([activity.loc[regulator], expression.loc[regulator]]).reshape(1, 2),
                              columns = ['Activity', 'Expression'], index = [regulator])

    sub_regulon = regulon.loc[regulator]
    up_reg = sub_regulon[(sub_regulon['MoA'] > 0)]
    down_reg = sub_regulon[(sub_regulon['MoA'] < 0)]

    empty_signature = pd.DataFrame(index = signature.index, columns = ['Up', 'Down']).fillna(0.0)
    empty_signature.loc[empty_signature.index.isin(up_reg.Target), 'Up'] = 1.0
    empty_signature.loc[empty_signature.index.isin(down_reg.Target), 'Down'] = 1.0

    up = empty_signature[(empty_signature['Up'] == 1)]['Up']
    down = empty_signature[(empty_signature['Down'] == 1)]['Down']

    pos_idx = [empty_signature.index.get_loc(x) for x in up.index.tolist()]
    neg_idx = [empty_signature.index.get_loc(x) for x in down.index.tolist()]

    if min(pos_idx) != 0:
        pos_idx.insert(0, 0)
    if min(neg_idx) != 0:
        neg_idx.insert(0, 0)
    if max(pos_idx) != len(signature):
        pos_idx.insert(-1, len(signature))
    if max(neg_idx) != len(signature):
        neg_idx.insert(-1, len(signature))

    fig = plt.figure(figsize = (14, 7))
    ax1 = fig.add_axes([0.1, 0.5, 0.8, 0.4], xticklabels = [], yticklabels = [])
    ax2 = fig.add_axes([0.1, 0.45, 0.8, 0.05], xticklabels = [], yticklabels = [])
    sub_ax1 = fig.add_axes([0.91, 0.8, 0.1, 0.1], xticklabels = [], yticklabels = [])
    sub_ax2 = fig.add_axes([0.91, 0.625, 0.1, 0.1], xticklabels = [], yticklabels = [])

    sns.distplot(pos_idx, rug = False, color = 'red', hist = False, kde_kws = {'cut': 0.0}, ax = ax1)

    sns.distplot(neg_idx, rug = False, color = 'blue', hist = False, kde_kws = {'cut': 0.0}, ax = ax1)

    sns.rugplot(pos_idx, color = 'r', **{'linewidth': .05, 'alpha': 0.25, 'height': 1}, ax = ax2)
    sns.rugplot(neg_idx, color = 'blue', **{'linewidth': .05, 'alpha': 0.25, 'height': 1}, ax = ax2)
    sns.heatmap(heat_frame, cmap = 'coolwarm', center = 0, ax = sub_ax1, annot = True, cbar = False, linecolor = 'k',
                linewidths = 1, yticklabels = False, annot_kws = {'size': 7})
    sub_ax1.set_ylabel('')

    sub_ax2.text(1, 7, 'p-value : {:.2E}'.format(Decimal(stats_frame.loc[regulator, 'pvalue'])), style = 'italic',
                 fontsize = 6)

    sub_ax2.text(1, 4, 'expression rank : {} / {}'.format(signature.index.get_loc(regulator) + 1, len(signature)),
                 style = 'italic', fontsize = 6)

    sub_ax2.text(1, 1, 'RES rank : {} / {}'.format(stats_frame.index.get_loc(regulator) + 1, stats_frame.shape[0]),
                 style = 'italic', fontsize = 6)
    sub_ax2.axis([0, 10, 0, 10])

    ax1.tick_params(bottom = False, left = True, direction = 'in')
    ax2.tick_params(bottom = False, left = False)
    ax2.axvline(x = signature.index.get_loc(regulator), color = 'k', linestyle = '--')

    ax2.text(0.1, -0.1, '(-) FC', style = 'italic', horizontalalignment = 'left', verticalalignment = 'top',
             fontsize = 6, fontweight = 'bold', transform = ax2.transAxes)
    ax2.text(0.9, -0.1, '(+) FC', style = 'italic', horizontalalignment = 'center', verticalalignment = 'top',
             fontsize = 6, fontweight = 'bold', transform = ax2.transAxes)

    ax2.text(0.5, -0.1, '0.0', style = 'italic', horizontalalignment = 'center', verticalalignment = 'top',
             fontsize = 6, fontweight = 'bold', transform = ax2.transAxes)

    ax2.set_xlabel('Ranked Expression Signature', style = 'italic')
    sub_ax1.tick_params(bottom = False, left = False, top = False)
    sub_ax2.tick_params(bottom = False, left = False, top = False)

    ax2.legend([Line2D([0], [0], color = 'white', visible = False),
                Line2D([0], [0], marker = '|', color = 'red', lw = 3.7, ls = 'None', markeredgewidth = 1.5),
                Line2D([0], [0], color = 'white', visible = False),
                Line2D([0], [0], marker = '|', color = 'blue', lw = 3.7, ls = 'None', markeredgewidth = 1.5),
                Line2D([0], [0], color = 'white', visible = False), Line2D([0], [0], color = 'white', lw = 3.7,
                                                                           path_effects = [pe.Stroke(linewidth = 4,
                                                                                                     foreground = 'k'),
                                                                                           pe.Normal()])],

               ["", "Activation (n = {})".format(len(up)), "", "Repression (n = {})".format(len(down)), "",
                "Background (n = {})".format(empty_signature.shape[0])], fontsize = 1, loc = 'center left',
               bbox_to_anchor = (1.0, 1.0), bbox_transform = ax2.transAxes, fancybox = False, shadow = False,
               title = 'Anticipatory Mode of Regulation', ncol = 1, prop = {'size': 7}, title_fontsize = 8,
               frameon = False)
    ax2.set_xticks([0, signature.shape[0] + 1])
    sns.despine(ax = sub_ax1)
    sns.despine(ax = sub_ax2, left = True, bottom = True, trim = True)
    sns.despine(ax = ax1, trim = True)
    sns.despine(ax = ax2, trim = True, top = False, left = True, bottom = False, right = True)

    ax1.set_title(regulator)
    plt.savefig('{}_res.pdf'.format(regulator), format = 'pdf', bbox_inches = 'tight')

priors_f = '/Users/estabroj/PycharmProjects/regulon_enrichment/data/causal-priors.txt'

def plot_label_distribution_gene_level_healthy_status(r, meta, mut, plot_dir, expr, total_enrichment, regulon):
    scores, pred_frame = split_fit_predictions(expr, total_enrichment, regulon, r, n_splits = 5, n_repeats = 10,
                          regressor = 'bayesridge')
    # cebpa_double = ['13-00342','13-00602','14-00034','14-00272','14-00777','15-00471','15-00786','16-00351','16-00880','16-01142','16-01219']

    rna_cebpa_double = ['13-00342', '13-00602', '14-00034', '15-00471', '15-00786']
    plt.close('all')
    plt.clf()
    infer_vals = pred_frame
    meta_filt = meta[meta.RNA_SeqID.isin(infer_vals.index)]
    meta_filt.set_index('RNA_SeqID', inplace = True)
    meta_filt = meta_filt.reindex(infer_vals.index)
    wt_clr = '0.29'
    control_clr = 'blue'
    mut_clrs = sns.light_palette('#C50000', reverse = True)

    fig, ax = plt.subplots(figsize = (7, 14))

    infer_means = infer_vals
    ens_gene = r
    filt_mut = mut[mut.symbol.isin([ens_gene])]

    ctd2_means = infer_vals.reindex(meta_filt[(meta_filt['PatientID'] != 'control')].index)

    control_means = infer_vals.reindex(meta_filt[(meta_filt['PatientID'] == 'control')].index)

    if np.all(infer_means >= 0):
        plt_ymin, plt_ymax = 0, max(np.max(infer_means.max()) * 1.09, 1)

    else:
        plt_ymax = np.max([np.max(np.absolute(infer_means.max())) * 1.09, 1.1])
        plt_ymin = -plt_ymax

    plt.ylim(plt_ymin, plt_ymax)
    plt_xmin, plt_xmax = plt.xlim()
    lbl_pad = (plt_ymax - plt_ymin) / 79

    mtype_stat = meta_filt.loc[ctd2_means.index].AML_Original_LabID.isin(filt_mut.seqid)
    kern_bw = (plt_ymax - plt_ymin) / 47

    for i in range(ctd2_means.shape[1]):
        ax = sns.kdeplot(ctd2_means.loc[~mtype_stat,str(i)], color = wt_clr, vertical = True, shade = True, alpha = 0.15,
                         linewidth = 0.1, bw = kern_bw, cut = 0, gridsize = 1000, label = 'Wild-Type')
        ax = sns.kdeplot(control_means.loc[:,str(i)], color = control_clr, vertical = True, shade = True, alpha = 0.15, linewidth = 0.1,
                         bw = kern_bw, gridsize = 1000, label = 'Control')

    ctd2_mut_means = list(zip(ctd2_means[mtype_stat].index, ctd2_means[mtype_stat].mean(axis=1)))
    ctd2_mut_std = dict(zip(ctd2_means[mtype_stat].index, ctd2_means[mtype_stat].std(axis = 1)))

    ctd2_cnt_means = list(zip(control_means.index, control_means.mean(axis=1),meta_filt.loc[control_means.index]['condition']))
    ctd2_cnt_std = dict(zip(control_means.index, control_means.std(axis = 1)))

    for i, (patient, val) in enumerate(ctd2_mut_means):
        if patient in rna_cebpa_double:
            print('plotting BAML CEBPA double mut patients')
            print(patient)
            plt_str = '{}'.format(patient)
            plt_clr = 'r'
            plt_lw = 1.7
            plt_lw_st = 3.0
            lw_st = ctd2_mut_std[patient]

            # TODO add linewidth scale and legend to increment std surrounding mutated samples
            # ax.axhline(y = val, xmin = 0, xmax = lw_st, c = 'white', alpha = 0.5, lw = plt_lw_st, path_effects=[pe.SimpleLineShadow(shadow_color='k'), pe.Normal()])
            ax.axhline(y = val, xmin = 0, xmax = lw_st, c = 'lightgray', alpha = 0.9, lw = plt_lw_st)
            ax.axhline(y = val, xmin = 0, xmax = plt_xmax * 0.22, c = plt_clr, ls = '--', lw = plt_lw)

            if i > 0 and ctd2_mut_means[i - 1][1] > (val - lbl_pad):
                txt_va = 'bottom'

            elif (i < (len(ctd2_mut_means) - 1) and ctd2_mut_means[i + 1][1] < (val + lbl_pad)):
                txt_va = 'top'

            else:
                txt_va = 'center'

            ax.text(plt_xmax * 0.32, val, r'$\bigstar$', size = 15, ha = 'left', color='gold',va = txt_va)


        else:
            print('plotting BAML patients')
            print(patient)
            plt_str = '{}'.format(patient)
            plt_clr = 'r'
            plt_lw = 1.7
            plt_lw_st = 3.0
            lw_st = ctd2_mut_std[patient]

            # TODO add linewidth scale and legend to increment std surrounding mutated samples
            # ax.axhline(y = val, xmin = 0, xmax = lw_st, c = 'white', alpha = 0.5, lw = plt_lw_st, path_effects=[pe.SimpleLineShadow(shadow_color='k'), pe.Normal()])
            ax.axhline(y = val, xmin = 0, xmax = lw_st, c = 'lightgray', alpha = 0.9, lw = plt_lw_st)
            ax.axhline(y = val, xmin = 0, xmax = plt_xmax * 0.22, c = plt_clr, ls = '--', lw = plt_lw)

        # if i > 0 and ctd2_mut_means[i - 1][1] > (val - lbl_pad):
        #     txt_va = 'bottom'
        #
        # elif (i < (len(ctd2_mut_means) - 1) and ctd2_mut_means[i + 1][1] < (val + lbl_pad)):
        #     txt_va = 'top'
        #
        # else:
        #     txt_va = 'center'
        #
        # ax.text(plt_xmax * 0.32, val, plt_str, size = 9, ha = 'left', va = txt_va)

    # cnt_color_map = {'cd34_pooled_technical' : 'darkorange', 'Individual_MNC':'darkgreen', 'Individual_CD34':'fuchsia'}
    cnt_color_map = {'cd34_pooled_technical': 'black', 'Individual_MNC': 'cyan', 'Individual_CD34': 'dimgrey'}

    for i, (patient, val, type) in enumerate(ctd2_cnt_means):
        print('plotting CTD2 patients')
        print(patient)
        plt_str = '{}'.format(patient)
        plt_clr = cnt_color_map[type]
        plt_lw = 1.7
        plt_lw_st = 3.0
        lw_st = ctd2_cnt_std[patient]

        # TODO add linewidth scale and legend to increment std surrounding mutated samples
        # ax.axhline(y = val, xmin = 0, xmax = lw_st, c = 'white', alpha = 0.5, lw = plt_lw_st, path_effects=[pe.SimpleLineShadow(shadow_color='k'), pe.Normal()])
        ax.axhline(y = val, xmin = 0, xmax = lw_st, c = 'lightgray', alpha = 0.9, lw = plt_lw_st)
        ax.axhline(y = val, xmin = 0, xmax = plt_xmax * 0.22, c = plt_clr, ls = '--', lw = plt_lw)


        # if i > 0 and ctd2_cnt_means[i - 1][1] > (val - lbl_pad):
        #     txt_va = 'bottom'
        #
        # elif (i < (len(ctd2_cnt_means) - 1) and ctd2_cnt_means[i + 1][1] < (val + lbl_pad)):
        #     txt_va = 'top'
        #
        # else:
        #     txt_va = 'center'
        #
        # ax.text(plt_xmax * 0.32, val, plt_str, size = 9, ha = 'left', va = txt_va)

    # calculate the accuracy of the mutation scores inferred across
    # validation runs in predicting mutation status
    # add annotation about the mutation scores' accuracy to the plot
    ax.text(ax.get_xlim()[1] * 0.91, plt_ymax * 0.82, size = 15, ha = 'right',
            s = "accuracy: {:2.3f} ".format(np.mean(scores)))

    # plt.xlabel('BAML', fontsize = 21, weight = 'semibold')

    plt.ylabel('Inferred {} Regulator Enrichment Score'.format(ens_gene), fontsize = 21, weight = 'semibold')

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

    print('Generating legend with tcga line')
    ax.legend([Patch(color = wt_clr, alpha = 0.36),
               Patch(color = control_clr, alpha = 0.36),

               Line2D([0], [0], color = 'white', visible = False),
               Line2D([0], [0], color = 'white', visible=False),
               Line2D([0], [0], color = mut_clrs[0], lw = 3.7, ls = '--'),
               Line2D([0], [0], color = cnt_color_map['cd34_pooled_technical'], lw = 3.7, ls = '--'),
               Line2D([0], [0], color = cnt_color_map['Individual_CD34'], lw = 3.7, ls = '--'),
               Line2D([0], [0], color = cnt_color_map['Individual_MNC'], lw = 3.7, ls = '--')],
              ["BAML Wild-Types",
               "BAML controls",
               '',
               '',
               "BAML {} mutants (n={})".format(ens_gene, mtype_stat.sum()),
               "CD34+ pooled reps (n={})".format(meta_filt.loc[control_means.index]['condition'].value_counts()['cd34_pooled_technical']),
               "Individual CD34+ (n={})".format(meta_filt.loc[control_means.index]['condition'].value_counts()['Individual_CD34']),
               "Individual MNC (n={})".format(meta_filt.loc[control_means.index]['condition'].value_counts()['Individual_MNC'])],
              fontsize = 13, loc = 'upper center', bbox_to_anchor = (0.5, -0.05), fancybox = True, shadow = True,
              ncol = 2)

    # scalebar = AnchoredSizeBar(ax.transData, .1, r'$BAML\ .1\ \sigma$', loc='lower left', color = 'whitesmoke', frameon = True,
    #                            size_vertical = .05, fontproperties = fontprops)
    # ax.add_artist(scalebar)


    ob = AnchoredHScaleBar(size = .1, label = r'$.1\ \sigma$', loc = 'upper left', frameon = False, pad = 0.0, sep = 4, color = "k")
    ax.add_artist(ob)
    fig.savefig(os.path.join(plot_dir, '{}_tfa_distribution.pdf'.format(ens_gene)), dpi = 300,
        bbox_inches = 'tight')

    plt.close()
    plt.gca()
    plt.clf()


def compare_moa(priors_f, regulator, regulon):
    plt.close('all')
    plt.clf()
    priors = pd.read_csv(priors_f, sep = '\t', index_col = 0, header = None)
    priors.columns = ['Type', 'Target', 'reference', 'site']
    sub_prior = priors.loc[regulator]
    prior_up_idx = sub_prior.Type.isin(['upregulates-expression'])
    prior_down_idx = sub_prior.Type.isin(['downregulates-expression'])
    sub_regul = regulon.loc[regulator]

    regul_up_idx = sub_regul['MoA'] > 0.0
    regul_down_idx = sub_regul['MoA'] < 0.0

    up_ = sub_prior.loc[prior_up_idx].Target.isin(sub_regul.loc[regul_up_idx].Target).sum()
    down_ = sub_prior.loc[prior_down_idx].Target.isin(sub_regul.loc[regul_down_idx].Target).sum()

    up_in = sub_prior.loc[prior_up_idx].Target.isin(sub_regul.loc[regul_down_idx].Target).sum()
    down_in = sub_prior.loc[prior_down_idx].Target.isin(sub_regul.loc[regul_up_idx].Target).sum()

    up_in_target = sub_prior.loc[prior_up_idx][
        sub_prior.loc[prior_up_idx].Target.isin(sub_regul.loc[regul_down_idx].Target)].Target
    down_in_target = sub_prior.loc[prior_down_idx][
        sub_prior.loc[prior_down_idx].Target.isin(sub_regul.loc[regul_up_idx].Target)].Target

    up_in_frame = sub_regul[(sub_regul.Target.isin(up_in_target))].loc[:, ['Target', 'MoA']]
    down_in_frame = sub_regul[(sub_regul.Target.isin(down_in_target))].loc[:, ['Target', 'MoA']]

    joined_inc = pd.concat([up_in_frame, down_in_frame])

    count_frame = pd.DataFrame.from_dict(
        {'Up': up_, 'Down': down_, 'Up-inconsistent': up_in, 'down-inconsistent': down_in}, orient = 'index',
        columns = ['counts'])
    count_frame.reset_index(inplace = True)
    count_frame.columns = ['consistency', 'Counts']

    plt.subplot(211)
    plt.title(regulator)
    ax1 = sns.barplot(x = 'consistency', y = 'Counts', data = count_frame, palette = ['r', 'b', 'b', 'r'])
    plt.subplot(212)
    clrs = ['red' if (x > 0) else 'blue' for x in joined_inc['MoA']]
    ax2 = sns.barplot(data = joined_inc, x = 'Target', y = 'MoA', palette = clrs)
    plt.ylim(1, -1)
    plt.xticks(rotation = 45, fontsize = 6)
    plt.tight_layout()
    plt.savefig('causal_priors_{}_check.pdf'.format(regulator), format = 'pdf', bbox_inches = 'tight')
    plt.close('all')
    plt.clf()


def generate_dag():
    dot = Digraph(comment = '{}'.format(regulator))
    regulators = regulon.index.unique().tolist()
    sub_regul = regulon.loc[regulator]
    sub_regul.loc[sub_regul.Target.isin(regulators)].Target


def generate_e_graph(meta, expr, total_enrichment, regulon, regulators):
    signature = expr.mean() / expr.std()

    control_samps = meta[(meta.condition.isin(['Individual_MNC']))]['RNA_SeqID']
    test_samps = meta[(meta.condition.isin(['ctd2']))]['RNA_SeqID']

    control = total_enrichment.loc[control_samps]
    test = total_enrichment.loc[test_samps]

    activity = test.mean() - control.mean()

    control_expr = expr.loc[control_samps]
    test_expr = expr.loc[test_samps]

    ttest_results_expr = st.ttest_ind(control_expr, test_expr)
    stats_frame_expr = pd.DataFrame(np.array([ttest_results_expr.pvalue, ttest_results_expr.statistic]),
                                    index = ['pvalue', 't-statistic'], columns = control_expr.columns).T

    weighted_ranks = (st.norm.ppf(stats_frame_expr['pvalue'] / 2.0, loc = 0, scale = 1)) * np.sign(
        stats_frame_expr['t-statistic'])
    weighted_ranks.sort_values(inplace = True)
    signature = signature.loc[weighted_ranks.index]

    expression = test_expr.mean() - control_expr.mean()

    ttest_results = st.ttest_ind(control, test)
    stats_frame = pd.DataFrame(np.array([ttest_results.pvalue, ttest_results.statistic]),
                               index = ['pvalue', 't-statistic'], columns = control.columns).T
    stats_frame.sort_values('pvalue', inplace = True)
    heat_frame = pd.DataFrame(np.array([activity.loc[regulators], expression.loc[regulators]]),
                              index = ['Activity', 'Expression'], columns = [regulators]).T

    cebpa = regulon.loc['CEBPA']

    cebpa_down = cebpa[(cebpa['MoA'] < 0)]
    cebpa_up = cebpa[(cebpa['MoA'] > 0)]

    e2f1 = regulon.loc['E2F1']

    e2f1_down = e2f1[(e2f1['MoA'] < 0)]
    e2f1_up = e2f1[(e2f1['MoA'] > 0)]

    mp = cebpa_down[cebpa_down.Target.isin(e2f1_up.Target)]
    pm = cebpa_up[cebpa_up.Target.isin(e2f1_down.Target)]

    empty_signature = pd.DataFrame(index = signature.index, columns = ['Up', 'Down']).fillna(0.0)
    empty_signature.loc[empty_signature.index.isin(mp.Target), 'Up'] = 1.0
    empty_signature.loc[empty_signature.index.isin(pm.Target), 'Down'] = 1.0

    up = empty_signature[(empty_signature['Up'] == 1)]['Up']
    down = empty_signature[(empty_signature['Down'] == 1)]['Down']

    pos_idx = [empty_signature.index.get_loc(x) for x in up.index.tolist()]
    neg_idx = [empty_signature.index.get_loc(x) for x in down.index.tolist()]

    if min(pos_idx) != 0:
        pos_idx.insert(0, 0)
    if min(neg_idx) != 0:
        neg_idx.insert(0, 0)
    if max(pos_idx) != len(signature):
        pos_idx.insert(-1, len(signature))
    if max(neg_idx) != len(signature):
        neg_idx.insert(-1, len(signature))

    fig = plt.figure()
    ax1 = fig.add_axes([0.1, 0.5, 0.8, 0.4], xticklabels = [], yticklabels = [])
    ax2 = fig.add_axes([0.1, 0.45, 0.8, 0.05], xticklabels = [], yticklabels = [])
    ax3 = fig.add_axes([0.1, 0.40, 0.8, 0.05], xticklabels = [], yticklabels = [])
    sub_ax1 = fig.add_axes([0.91, 0.8, 0.1, 0.1])
    sub_ax2 = fig.add_axes([0.1, 0.9, 0.3, 0.1], xticklabels = [], yticklabels = [])

    sns.distplot(pos_idx, rug = False, color = 'red', hist = False, kde_kws = {'cut': 0.0}, ax = ax1)

    sns.distplot(neg_idx, rug = False, color = 'blue', hist = False, kde_kws = {'cut': 0.0}, ax = ax1)

    sns.rugplot(pos_idx, color = 'r', **{'linewidth': .00001, 'alpha': 0.05, 'height': 1}, ax = ax2)
    sns.rugplot(neg_idx, color = 'blue', **{'linewidth': .00001, 'alpha': 0.05, 'height': 1}, ax = ax3)
    sns.heatmap(heat_frame, cmap = 'coolwarm', center = 0, ax = sub_ax1, annot = True, cbar = False, linecolor = 'k',
                linewidths = 1, cbar_kws = {'label': ''}, annot_kws = {'size': 7})
    sub_ax1.set_ylabel('')
    sub_ax2.text(1, 7, 'p-value : {:.2E}'.format(Decimal(stats_frame.loc[regulator, 'pvalue'])), style = 'italic',
                 fontsize = 6, fontweight = 'bold')

    sub_ax2.text(1, 4, 'expression rank : {} / {}'.format(signature.index.get_loc(regulator) + 1, len(signature)),
                 style = 'italic', fontsize = 6, fontweight = 'bold')

    sub_ax2.text(1, 1, 'RES rank : {} / {}'.format(stats_frame.index.get_loc(regulator) + 1, stats_frame.shape[0]),
                 style = 'italic', fontsize = 6, fontweight = 'bold')
    sub_ax2.axis([0, 10, 0, 10])

    ax1.tick_params(bottom = False, left = False)
    ax2.tick_params(bottom = False, left = False)
    ax2.axvline(x = signature.index.get_loc(regulator), color = '#c46a21', linestyle = '--')
    ax3.tick_params(bottom = False, left = False)
    ax3.axvline(x = signature.index.get_loc(regulator), color = '#c46a21', linestyle = '--')

    ax2.text(0.99, 0.5, r'{} $CEBPA- \bigcap\  E2F1+$'.format(len(up)), style = 'italic', horizontalalignment = 'right',
             verticalalignment = 'center', fontsize = 4, transform = ax2.transAxes)
    ax3.text(0.1, -0.1, '(-) FC', style = 'italic', horizontalalignment = 'left', verticalalignment = 'top',
             fontsize = 6, fontweight = 'bold', transform = ax3.transAxes)
    ax3.text(0.9, -0.1, '(+) FC', style = 'italic', horizontalalignment = 'center', verticalalignment = 'top',
             fontsize = 6, fontweight = 'bold', transform = ax3.transAxes)
    ax3.text(0.99, 0.5, r'{} $CEBPA+ \bigcap\  E2F1-$'.format(len(down)), style = 'italic',
             horizontalalignment = 'right', verticalalignment = 'center', fontsize = 4, transform = ax3.transAxes)
    ax3.text(0.5, -0.1, '0.0', style = 'italic', horizontalalignment = 'center', verticalalignment = 'top',
             fontsize = 6, fontweight = 'bold', transform = ax3.transAxes)

    ax3.set_xlabel('Ranked Expression Signature', style = 'italic')
    sub_ax1.tick_params(bottom = False, top = False, right = False)
    sub_ax2.tick_params(bottom = False, left = False, top = False)

    ax2.legend([Line2D([0], [0], color = 'white', visible = False),
                Line2D([0], [0], marker = '|', color = 'red', lw = 3.7, ls = 'None', markeredgewidth = 1.5),
                Line2D([0], [0], color = 'white', visible = False),
                Line2D([0], [0], marker = '|', color = 'blue', lw = 3.7, ls = 'None', markeredgewidth = 1.5),
                Line2D([0], [0], color = 'white', visible = False), Line2D([0], [0], color = 'white', lw = 3.7,
                                                                           path_effects = [pe.Stroke(linewidth = 4,
                                                                                                     foreground = 'k'),
                                                                                           pe.Normal()])],

               ["", "Activation", "", "Repression", "", "Background", ], fontsize = 2, loc = 'center left',
               bbox_to_anchor = (1.0, 0.25), bbox_transform = ax2.transAxes, fancybox = False, shadow = False,
               title = 'Anticipatory Mode of Regulation', ncol = 1, prop = {'size': 3}, title_fontsize = 4,
               frameon = False)

    ## anticipatory regulation
    ## ETF1 targets negatively regulated

    plt.title(regulator)
    plt.savefig('{}_{}_res.pdf'.format(regulators[0], regulators[1]), format = 'pdf', bbox_inches = 'tight')


def find_recursive_shadow(regulon):
    shadow_dict = {}
    regulators = regulon.index.unique().tolist()
    for i in tqdm(range(len(regulators))):
        for j in tqdm(range(i + 1, len(regulators))):
            sub_regul_i = regulon.loc[regulators[i]].set_index('Target')
            sub_regul_j = regulon.loc[regulators[j]].set_index('Target')
            reg_i = regulators[i]
            reg_j = regulators[j]
            comparison = '{}_{}'.format(reg_i, reg_j)
            percent = 100 * (len(set(sub_regul_i.index).intersection(sub_regul_j.index)) / len(
                set(sub_regul_i.index).union(sub_regul_j.index)))
            shared_genes = set(sub_regul_i.index).intersection(sub_regul_j.index)
            rho = st.spearmanr(sub_regul_i.loc[shared_genes]['MoA'], sub_regul_j.loc[shared_genes]['MoA'])[0]
            if comparison not in shadow_dict:
                shadow_dict[comparison] = [reg_i, reg_j, percent, rho]


def load_ctd2():
    path = '/Users/estabroj/PycharmProjects/regulon_enrichment/experiments/CTD2/data/misc'
    regulon = read_pickle(os.path.join(path, 'CTD2_bolstered_regulon.pkl'))
    expr = read_pickle(os.path.join(path, 'CTD2_robust_True_False_0.4_True_frame.pkl'))
    shadow = pd.read_csv(os.path.join(path, 'ctd2_shadow_regulon.txt'), sep = '\t', index_col = 0)
    filt_shadow = shadow[(shadow['percent'] > 50)]
    filt_shadow.columns = ['primary_regulator', 'secondary_regulator', 'percent', 'correlation']
    cohort = 'CTD2'
    meta = pd.read_csv(os.path.join(path, 'ctd2_meta_with_healthy_status.txt'), sep = '\t', index_col = 0)
    total_enrichment = read_pickle(os.path.join(path, 'CTD2_total_enrichment.pkl'))
    stats_frame = pd.read_csv(os.path.join(path, 'CTD2_stats_frame.txt'), sep = '\t', index_col = 0)

    return total_enrichment, meta, cohort, shadow, filt_shadow, expr, regulon, stats_frame


def find_shadow_pairs(filt_shadow, regulator, regulators, rho = .70):
    shadow_pair = filt_shadow[
        (filt_shadow['primary_regulator'] == regulator) | (filt_shadow['secondary_regulator'] == regulator)]
    shadow_pair_sig = shadow_pair[
        (shadow_pair['primary_regulator'].isin(regulators)) & (shadow_pair['secondary_regulator'].isin(regulators))]
    shadow_pair_rho = shadow_pair_sig[(abs(shadow_pair_sig['correlation']) >= rho)]
    u_set = set(shadow_pair_rho.primary_regulator).union(shadow_pair_rho.secondary_regulator)
    u_set.remove(regulator)

    return u_set


def find_shadow_pairs_recursive(filt_shadow, regulator, regulators, rho = .70):
    idxs = set()
    shadow_pair = filt_shadow[
        (filt_shadow['primary_regulator'] == regulator) | (filt_shadow['secondary_regulator'] == regulator)]
    shadow_pair_sig = shadow_pair[
        (shadow_pair['primary_regulator'].isin(regulators)) & (shadow_pair['secondary_regulator'].isin(regulators))]

    shadow_pair_rho = shadow_pair_sig[(abs(shadow_pair_sig['correlation']) >= rho)]
    idxs.update(shadow_pair_rho.index.tolist())
    u_set = set(shadow_pair_rho.primary_regulator).union(shadow_pair_rho.secondary_regulator)
    u_set.remove(regulator)

    for r in u_set:
        print(r)
        print(u_set)
        shadow_pair = filt_shadow[(filt_shadow['primary_regulator'] == r) | (filt_shadow['secondary_regulator'] == r)]
        shadow_pair_sig = shadow_pair[
            (shadow_pair['primary_regulator'].isin(regulators)) & (shadow_pair['secondary_regulator'].isin(regulators))]

        shadow_pair_rho = shadow_pair_sig[(abs(shadow_pair_sig['correlation']) >= rho)]
        idxs.update(shadow_pair_rho.index.tolist())

    return idxs


def identify_coregulation(expr, regulon, regulator, filt_shadow, cohort, stats_frame):
    quant_nes = load_quantile(regulon, expr, cohort)
    regulators = stats_frame[(stats_frame.pvalue <= 0.05)].index
    s_pairs = find_shadow_pairs(filt_shadow, regulator, regulators, rho = .1)
    regulator_regulon = regulon.loc[regulator]
    control_samps = meta[(meta.condition.isin(['Individual_MNC']))]['RNA_SeqID']
    test_samps = meta[(meta.condition.isin(['ctd2']))]['RNA_SeqID']

    shadow_dict = {}

    for sg in s_pairs:
        shadow_regulon = regulon.loc[sg]
        unique_shadow = shadow_regulon[~(shadow_regulon.Target.isin(regulator_regulon.Target))]
        unique_regulator = regulator_regulon[~(regulator_regulon.Target.isin(shadow_regulon.Target))]
        joined_regulon = pd.concat([unique_shadow, unique_regulator])
        if sum(joined_regulon.index.value_counts() > 1) == 2:
            regulator_list = [regulator, sg]

            print('--- Calculating regulon enrichment scores ---')
            nes_list = list(
                map(functools.partial(score_enrichment, expr = expr, regulon = joined_regulon, quant_nes = quant_nes),
                    tqdm(regulator_list)))
            shadow_enrichment = pd.concat(nes_list, axis = 1)

            control = shadow_enrichment.loc[control_samps]
            test = shadow_enrichment.loc[test_samps]
            ttest_results = st.ttest_ind(control, test)
            shadow_stats_frame = pd.DataFrame(np.array([ttest_results.pvalue, ttest_results.statistic]),
                                              index = ['pvalue', 't-statistic'], columns = control.columns).T
            shadow_stats_frame.sort_values('pvalue', inplace = True)

            l10 = np.log10(shadow_stats_frame['pvalue'])
            pde = l10[1] - l10[0]

            bool_pvalues = stats_frame.loc[regulator_list]['pvalue'] < shadow_stats_frame.loc[regulator_list]['pvalue']

            if sg not in shadow_dict:
                shadow_dict[sg] = [bool_pvalues]
        else:
            continue

    return shadow_dict


def bin_edge_types(shadow_pair_dict):
    # increase edge - directed edge pointing from node A to node B
    # - >
    increase_edge = []

    # decrease edge - directed edge pointing to node A from node B
    # ] -
    decrease_edge = []

    # synergy edge - directed edge pointing from node A to/from node B
    # < - >
    synergy_edge = []

    # null edge - no directed edge pointing from node A to/from node B
    null_edge = []

    for k in shadow_pair_dict:
        if shadow_pair_dict[k][0].values.tolist() == [True, False]:
            decrease_edge.append(tuple(shadow_pair_dict[k][0].index.tolist())[::-1])
        if shadow_pair_dict[k][0].values.tolist() == [False, True]:
            increase_edge.append(tuple(shadow_pair_dict[k][0].index.tolist()))
        if shadow_pair_dict[k][0].values.tolist() == [True, True]:
            synergy_edge.append(tuple(shadow_pair_dict[k][0].index.tolist()))
        if shadow_pair_dict[k][0].values.tolist() == [False, False]:
            null_edge.append(tuple(shadow_pair_dict[k][0].index.tolist()))
    return increase_edge, decrease_edge, synergy_edge, null_edge


def identify_coregulation_idxs(expr, regulon, filt_shadow, cohort, stats_frame, idxs):
    quant_nes = load_quantile(regulon, expr, cohort)
    s_pairs = filt_shadow.loc[idxs, ['primary_regulator', 'secondary_regulator']].values
    control_samps = meta[(meta.condition.isin(['Individual_MNC']))]['RNA_SeqID']
    test_samps = meta[(meta.condition.isin(['ctd2']))]['RNA_SeqID']
    shadow_dict = {}

    for paired in s_pairs:
        regulator = paired[0]
        sg = paired[1]
        edge = '{}_{}'.format(regulator, sg)
        regulator_regulon = regulon.loc[regulator]
        shadow_regulon = regulon.loc[sg]
        unique_shadow = shadow_regulon[~(shadow_regulon.Target.isin(regulator_regulon.Target))]
        unique_regulator = regulator_regulon[~(regulator_regulon.Target.isin(shadow_regulon.Target))]
        joined_regulon = pd.concat([unique_shadow, unique_regulator])
        regulator_list = [regulator, sg]

        print('--- Calculating regulon enrichment scores ---')
        nes_list = list(
            map(functools.partial(score_enrichment, expr = expr, regulon = joined_regulon, quant_nes = quant_nes),
                tqdm(regulator_list)))
        shadow_enrichment = pd.concat(nes_list, axis = 1)

        control = shadow_enrichment.loc[control_samps]
        test = shadow_enrichment.loc[test_samps]
        ttest_results = st.ttest_ind(control, test)
        shadow_stats_frame = pd.DataFrame(np.array([ttest_results.pvalue, ttest_results.statistic]),
                                          index = ['pvalue', 't-statistic'], columns = control.columns).T
        shadow_stats_frame.sort_values('pvalue', inplace = True)

        bool_pvalues = stats_frame.loc[regulator_list]['pvalue'] < shadow_stats_frame.loc[regulator_list]['pvalue']

        if edge not in shadow_dict:
            shadow_dict[edge] = [bool_pvalues]

    return shadow_dict


def simple_graph():
    import networkx as nx
    import matplotlib.pyplot as plt

    total_graph = nx.from_pandas_edgelist(filt_shadow, 'primary_regulator', 'secondary_regulator',
                                          edge_attr = 'correlation', create_using = nx.DiGraph())
    graph = total_graph.subgraph(nx.shortest_path(total_graph.to_undirected(), 'CEBPA'))
    pos = nx.nx_agraph.graphviz_layout(graph, prog = 'sfdp', args = '-Goverlap=false')

    regulators = set(filt_shadow.primary_regulator) | set(filt_shadow.secondary_regulator)

    val_map = test.loc[:, list(regulators)].median().to_dict()

    values = [val_map.get(node, 0.0) for node in graph.nodes()]

    labels = {node: node for node in graph.nodes()}
    nx.set_node_attributes(graph, 'labels', labels.values())
    nodes = nx.draw_networkx_nodes(graph, pos, cmap = plt.get_cmap('plasma'), node_color = values, node_size = 400,
                                   alpha = 1, with_labels = True)
    nx.draw_networkx_labels(graph, pos, labels, font_size = 7)
    edge_idx = ['{}_{}'.format(u, v) for u, v in graph.edges()]
    colors = filt_shadow.loc[edge_idx]['correlation'].values.tolist()
    edges = nx.drawing.nx_pylab.draw_networkx_edges(graph, pos, edge_cimport pandas as pd

import os
import subprocess
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats

DATA_PATH = os.getcwd()


def load_phospho_data_raw_plex(xlsx_phospho_data='CEDAR_phospho_peptide_abundances_sent.xlsx', sheet_name_phos = 'Sheet1'):
    """ Load beta4 phospho data

    Args:
        xlsx_phospho_data (str): title of xlsx spreadsheet of phospho dataset
        sheet_name_phos (str): xlsx sheet to read into pandas frame
    Returns:
        phospho_data (:obj: `pandas DataFrame`) : pandas DataFrame containing processed phospho data
    """
    data = pd.read_excel(os.path.join(DATA_PATH, xlsx_phospho_data), sheet_name = sheet_name_phos)

    return data


def split_by_plex(phospho_raw):
    phospho_raw.set_index(['Peptide', 'FirstOfProtein'], inplace = True)
    cols = ['KC_C1891', 'KC_B4_C1844', 'KC_B4_C1907', 'WT_C1893', 'KC_C1818', 'KC_C1903', 'KC_B4_C1889',
            'WT_C1908_rep1', 'KC_C1898', 'KC_C1848', 'KC_B4_C1935', 'WT_C1909', 'KC_C1793', 'KC_B4_C1882',
            'KC_B4_C1847', 'WT_C1908_rep2']
    ref1 = phospho_raw.iloc[:, :9]
    ref1 = ref1.subtract(ref1.median())
    ref1_ = ref1.Ref1.to_frame()
    ref1 = ref1.iloc[:, 1:]
    ref1_normed = (ref1 - ref1_.values)


    ref2 = phospho_raw.iloc[:, 9:]
    ref2 = ref2.subtract(ref2.median())
    ref2_ = ref2.Ref2.to_frame()
    ref2 = ref2.iloc[:, 1:]
    ref2_normed = (ref2 - ref2_.values)

    joined_phospho = pd.concat([ref1_normed, ref2_normed],axis=1)



def load_phospho_data_raw(xlsx_phospho_data='CEDAR_phospho.xlsx', sheet_name_phos = 'redundant_site'):
    """ Load beta4 phospho data

    Args:
        xlsx_phospho_data (str): title of xlsx spreadsheet of phospho dataset
        sheet_name_phos (str): xlsx sheet to read into pandas frame
    Returns:
        phospho_data (:obj: `pandas DataFrame`) : pandas DataFrame containing processed phospho data
    """
    data = pd.read_excel(os.path.join(DATA_PATH, xlsx_phospho_data), sheet_name = sheet_name_phos)

    return data


def load_total_protein_raw(xlsx_prot_data='CEDAR_Global.xlsx', sheet_name_prot = 'Global_ratio'):
    """ Load beta4 proteomics data

    Args:
        xlsx_prot_data (str): title of xlsx spreadsheet of proteomics dataset
        sheet_name_prot (str): xlsx sheet to read into pandas frame
    Returns:
        prot_data (:obj: `pandas DataFrame`) : pandas DataFrame containing processed proteomics data

    """


    data = pd.read_excel(os.path.join(DATA_PATH, xlsx_prot_data), sheet_name = sheet_name_prot)

    return data


def load_phospho_data(xlsx_phospho_data='CEDAR_phospho.xlsx', sheet_name_phos = 'redundant_site'):
    """ Load beta4 phospho data

    Args:
        xlsx_phospho_data (str): title of xlsx spreadsheet of phospho dataset
        sheet_name_phos (str): xlsx sheet to read into pandas frame
    Returns:
        phospho_data (:obj: `pandas DataFrame`) : pandas DataFrame containing processed phospho data
    """

    data = pd.read_excel(os.path.join(DATA_PATH, xlsx_phospho_data), sheet_name = sheet_name_phos)
    data = data.loc[~data.Protein_Residue_Num.isna()]
    na_genes = data.loc[data['Gene names'].isna(), 'Protein'].str.split('_', expand = True)[0]
    data.loc[data['Gene names'].isna(), 'Gene names'] = na_genes
    data['Gene names'].fillna('unknown', inplace = True)

    # s = data['Gene names'].str.split(" ").apply(pd.Series, 1).stack()
    # s.index = s.index.droplevel(-1)
    # s.name = 'Symbols'
    # del data['Gene names']
    # joined_data = data.join(s)

    s = data['Gene names'].str.split(" ",expand=True)[0].str.upper()
    s.name = 'Symbols'
    joined_data = data.join(s)

    joined_data['Position'] = joined_data['Protein_Residue_Num'].astype(int)
    joined_data['Sites'] = joined_data['Residue'] + joined_data['Position'].astype(str)
    joined_data['ID'] = joined_data['Symbols'].str.upper() + '_' + joined_data['Sites']
    joined_data['Effect'] = ''

    drop_dups = joined_data.loc[joined_data.notnull().sum(1).groupby(joined_data.ID).idxmax()]
    drop_dups.sort_values('ID', inplace = True)
    phospho_data = drop_dups.drop_duplicates()
    phospho_data.reset_index(drop = True, inplace = True)
    phospho_data['Symbols'] = phospho_data.Symbols.str.upper()

    s = phospho_data.Uniprot.str.split("; ").apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'Leading Uniprot ID'
    phospho_data = phospho_data.join(s)

    return phospho_data


def load_total_protein(xlsx_prot_data='CEDAR_Global.xlsx', sheet_name_prot = 'Global_ratio'):
    """ Load beta4 proteomics data

    Args:
        xlsx_prot_data (str): title of xlsx spreadsheet of proteomics dataset
        sheet_name_prot (str): xlsx sheet to read into pandas frame
    Returns:
        prot_data (:obj: `pandas DataFrame`) : pandas DataFrame containing processed proteomics data

    """

    ortho_genes = load_ortho_map()
    data = pd.read_excel(os.path.join(DATA_PATH, xlsx_prot_data), sheet_name = sheet_name_prot)
    data = data.loc[~data['gene name'].isna()]

    data['gene name'].fillna('unknown', inplace = True)

    ## duplicate values need to be dropped

    # s = data['gene name'].str.split(" ").apply(pd.Series, 1).stack()
    # s.index = s.index.droplevel(-1)
    # s.name = 'Symbols'
    # del data['gene name']

    s = data['gene name'].str.split(" ",expand=True)[0].str.upper()
    s.name = 'Symbols'
    joined_data = data.join(s)

    drop_dups = joined_data.loc[joined_data.notnull().sum(1).groupby(joined_data['Symbols']).idxmax()]
    drop_dups.sort_values('Symbols', inplace = True)
    prot_data = drop_dups.drop_duplicates()
    prot_data.reset_index(drop = True, inplace = True)
    prot_data['ID'] = prot_data['Symbols']
    prot_data['Sites'] = ''
    prot_data['Effect'] = ''
    prot_data['Symbols'] = prot_data.Symbols.str.upper()
    prot_data = prot_data[(prot_data['Symbols'].isin(ortho_genes))]
    prot_data_cols = sorted(prot_data.columns[6:-4])
    base_lab = ['Symbols','ID','Sites','Effect']
    joined_prot_labels = base_lab + prot_data_cols
    ordered_prot = prot_data.loc[:,joined_prot_labels]

    return ordered_prot


def load_ortho_map(ortho_data='mart_export.txt'):
    """ Load Ensembl gene list of known orthologue mapping between Mouse and Human:
        http://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/

    Args:
        ortho_data (str): file exported Ensembl with known mouse-human mapping

    Returns:
       ortho_genes (numpy.ndarray): array of conserved genes between mouse and human

    """

    data = pd.read_csv(os.path.join(DATA_PATH, ortho_data),sep='\t')
    filt_data = data[(data['Mouse orthology confidence [0 low, 1 high]'] == 1)]
    ortho_genes = filt_data['Mouse gene name'].str.upper().unique()

    return ortho_genes


def load_mouse_mappings(map_file='Mouse_to_Human_mapping.txt.gz'):
    """ Loads Dapple2 conserved phosphosite mapping from mouse to human

    Args:
        map_file (str): File handle of mouse to human phosphosite mapping

    Returns:
        map_frame (:obj: `pandas DataFrame`): UniProt and gene symbol mapping between mouse and human

    """

    map_frame = pd.read_csv(os.path.join(DATA_PATH, map_file),sep='\t',index_col=0)

    return map_frame


def load_meta(meta_f='phospho_beta_meta.txt'):
    """ Generates metadata file based on the unique sample headers in the phospho-proteomics dataset

    Args:
        phospho_data (:obj: `pandas DataFrame`) : pandas DataFrame containing processed phospho data
        prot_data (:obj: `pandas DataFrame`) : pandas DataFrame containing processed proteomics data

    Returns:
        merged_meta (:obj: `pandas DataFrame`) : frame containing metadata associated with both datasets:
            phosphoID, Time, Phase, protID

    """

    meta = pd.read_csv(os.path.join(DATA_PATH, meta_f),sep='\t',index_col=0)

    return meta


def merge_data(phospho_data, prot_data, map_frame):
    """ Merge phosphorylation and proteomics dataset. Integrate fuzzy site matching in dataset accounting for sites that
        do not map perfectly, but show to be conserved and the difference between their site number and the conserved
        site is 1. i.e. identified site: S141 & conserved site: S140. S140 will be included in the dataset and will
        represent the values identified for S141.

    Args:
        phospho_data (:obj: `pandas DataFrame`) : pandas DataFrame containing processed phospho data
        prot_data (:obj: `pandas DataFrame`) : pandas DataFrame containing processed proteomics data
        map_frame (:obj: `pandas DataFrame`): UniProt and gene symbol mapping between mouse and human

    Returns:
        phospho_prot (:obj: `pandas DataFrame`) : pandas DataFrame containing processed phospho-proteomics data

    """
    map_frame.reset_index(drop = True,inplace=True)
    map_frame['ID'] = map_frame['Gene'] + '_' + map_frame['Hit site']
    map_frame['prot_ID'] = map_frame['Hit accession'] + '_' + map_frame['Hit site']

    phospho_data['ID'] = phospho_data['ID'].str.upper()
    phospho_data['prot_ID'] = phospho_data['Leading Uniprot ID'] + '_' + phospho_data['Sites']

    matched_phospho = perfect_match(phospho_data, map_frame)

    matched_ids = matched_phospho.columns[matched_phospho.columns.str.contains('_')].tolist()
    matched_ids.pop(matched_ids.index('Protein_Residue_Num'))
    matched_ids.pop(matched_ids.index('prot_ID'))
    matched_ids = sorted(matched_ids)

    phospho_filt = ['ID','Symbols','Query site','Effect'] + matched_ids
    phospho_filt_data = matched_phospho.loc[:,phospho_filt].drop_duplicates()
    phospho_filt_data.rename(index = str, columns = {"Query site": "Sites"},inplace=True)

    fuzz_frame = fuzzy_match(phospho_data, map_frame)
    fuzz_filt = ['ID','Symbols','Query site','Effect'] + matched_ids
    fuzz_filt_data = fuzz_frame.loc[:,fuzz_filt]
    fuzz_filt_data.rename(index = str, columns = {"Query site": "Sites"},inplace=True)

    prot_map = {'WT_C1908_rep1':'plex1_WT_C1908','WT_C1908_rep2':'plex2_WT_C1908'}
    matched_prot_ids = list(map(lambda x: x if x not in prot_map else prot_map[x], matched_ids))
    prot_filt = ['ID','Symbols','Sites','Effect'] + matched_prot_ids
    rename = ['ID','Symbols','Sites','Effect'] + matched_ids
    prot_filt_data = prot_data.loc[:,prot_filt]
    prot_filt_data.columns = rename

    phospho_prot = pd.concat([phospho_filt_data, fuzz_filt_data, prot_filt_data], sort=False).sort_values(by='Symbols').drop_duplicates()

    return phospho_prot


def perfect_match(phospho_data, map_frame):
    """ Identify perfect match identified sites in phospho dataset and conserved Dapple2 mapping

    Args:
        phospho_data (:obj: `pandas DataFrame`) : pandas DataFrame containing processed phospho data
        map_frame (:obj: `pandas DataFrame`): UniProt and gene symbol mapping between mouse and human

    Returns:
        joined_match (:obj: `pandas DataFrame`) : pandas DataFrame containing processed phospho data with conserved
            phosphosites from mouse to human

    """
    phos_match = phospho_data[phospho_data.ID.isin(map_frame.ID)]
    map_match = map_frame.loc[map_frame.ID.isin(phos_match['ID'])]
    dup_idx = map_match.loc[:,['Hit site','RBH?','ID']].duplicated()
    map_frame_match = map_match.loc[~dup_idx].sort_values('RBH?',ascending=False)
    de_dup_idx = map_frame_match.loc[:, ['Hit site','ID']].duplicated()
    map_frame_de_dup = map_frame_match.loc[~de_dup_idx]

    joined_match = phos_match.merge(map_frame_de_dup.loc[:, ['ID', 'RBH?', 'Query site', 'Hit site']], on = 'ID')

    return joined_match


def fuzzy_match(phospho_data, map_frame):
    """ Identify fuzzy match phosphorylation sites in phospho dataset and conserved Dapple2 mapping

    Args:
        phospho_data (:obj: `pandas DataFrame`) : pandas DataFrame containing processed phospho data
        map_frame (:obj: `pandas DataFrame`): UniProt and gene symbol mapping between mouse and human

    Returns:
        joined_match (:obj: `pandas DataFrame`) : pandas DataFrame containing processed phospho data with conserved
            phosphosites from mouse to human

    """
    idx_p = phospho_data.prot_ID.isin(map_frame.prot_ID)
    idx_g = phospho_data.ID.isin(map_frame.ID)

    alt_idx = (idx_g + idx_p)
    filt_map = map_frame[map_frame.Gene.isin(phospho_data.loc[~alt_idx]['Symbols'].unique())]
    filt_phos = phospho_data.loc[~alt_idx]

    map_dict = {gene: group for gene,group in filt_map.groupby('Gene')}
    phos_dict = {gene:group for gene,group in filt_phos.groupby('Symbols')}

    fuzz_list = []
    genes = list(set(phos_dict.keys()).intersection(set(map_dict.keys())))
    for gene in genes:
        phos_filt = phos_dict[gene]
        map_filt = map_dict[gene].loc[map_dict[gene]['Hit site'].str.contains('|'.join(phos_dict[gene]['Sites'].str.slice(0, -1).values))]
        map_filt['Amino acid'] = map_filt['Hit site'].str.slice(0,1)
        map_filt['Position'] = map_filt['Hit site'].str.slice(1).astype(int)

        mapped = []
        for pindex, prow in phos_filt.iterrows():
            for mindex, mrow in map_filt.iterrows():
                if mrow['Amino acid'] == prow['Residue']:
                    if abs(mrow['Position'] - prow['Position']) <= 1:
                        prow['Query site'] = mrow['Query site']
                        prow['Hit site'] = mrow['Hit site']
                        prow['RBH?'] = mrow['RBH?']
                        mapped.append(prow.to_frame().T)
        if len(mapped) > 0:
            if len(mapped) == 1:
                fuzz_list.append(mapped[0])
            else:
                sub_fuzz_frame = pd.concat(mapped)
                fuzz_list.append(sub_fuzz_frame)
        else:
            continue

    fuzz_frame = pd.concat(fuzz_list)

    return fuzz_frame


def generate_parameter_file(relnm, test_samps, control_samps, value_transformation = 'significant-change-of-mean',
                            fdr_threshold = '0.1', site_match = '5', site_effect = '5'):
    """ Generate the required CausalPath input parameters file that associates the ProteomicsData file with the
        parameters for the analysis

    Args:
        relnm (str): out path for analysis files
        test_samps (list): list of samples to be used as the test group or contrast in the analysis
        control_samps (list): list of samples to be used as the control group or baseline in the analysis
        value_transformation (str): means of detecting changes between the test and control group
        fdr_threshold (str): False discovery rate threshold
        site_match (str): Associate sites with known site if site falls within n=site_match of known site
        site_effect (str): Associate sites with known activity of site if site falls within n=site_match of known site

    """
    initial_lines = ["proteomics-values-file = ProteomicsData.txt", "id-column = ID", "symbols-column = Symbols",
                     "sites-column = Sites", "effect-column = Effect", "do-log-transform = false"]
    out_f = open(os.path.join(relnm, 'parameters.txt'), 'w')
    for v in initial_lines:
        out_f.write(v + '\n')
    out_f.write('fdr-threshold-for-data-significance = ' + fdr_threshold + ' protein' + '\n')
    out_f.write('fdr-threshold-for-data-significance = ' + fdr_threshold + ' phosphoprotein' + '\n')
    out_f.write('value-transformation = ' + value_transformation + '\n')
    out_f.write('minimum-sample-size = 3' + '\n')
    out_f.write('calculate-network-significance = true' + '\n')
    out_f.write('pool-proteomics-for-fdr-adjustment = false' + '\n')
    out_f.write('permutations-for-significance = 1000' + '\n')
    out_f.write('color-saturation-value = 2' + '\n')
    out_f.write('site-match-proximity-threshold = ' + site_match + '\n')
    out_f.write('site-effect-proximity-threshold = ' + site_effect + '\n')
    out_f.write('show-insignificant-data = false'+ '\n')

    if value_transformation == 'correlation':
        for c in control_samps:
            out_f.write('value-column = ' + c + '\n')
    else:

        for c in control_samps:
            out_f.write('control-value-column = ' + c + '\n')
        for t in test_samps:
            out_f.write('test-value-column = ' + t + '\n')

    out_f.close()


def generate_proteomics_data(sub_data, relnm):
    """ Write out phospho-proteomices datafile to designated output directory

    Args:
        sub_data (:obj: `pandas DataFrame`) : pandas DataFrame containing processed phospho data with conserved
            phosphosites from mouse to human
        relnm (str): out path for analysis files

    """

    sub_data.iloc[:, 4:] = sub_data.iloc[:, 4:].fillna('NaN')
    sub_data.to_csv(path_or_buf=os.path.join(relnm, 'ProteomicsData.txt'),sep='\t')


def generate_data_files(phospho_prot, merged_meta, **kwargs):
    """

    Args:
        phospho_prot:
        merged_meta:
        **kwargs:

    Returns:

    """
    print(kwargs)
    print('Generating Time dependent ratio contrast')
    subset_meta = merged_meta[merged_meta['Condition'].isin(kwargs['Condition'])]
    samps_idx = ['ID','Symbols','Sites','Effect'] + subset_meta.index.tolist()
    pre_rx = subset_meta[(subset_meta['Condition'] == kwargs['Condition'][0])].index
    post_rx = subset_meta[(subset_meta['Condition'] == kwargs['Condition'][1])].index
    sub_data = phospho_prot.loc[:, samps_idx]
    sub_data.set_index('ID',inplace=True)

    return sub_data, pre_rx, post_rx


def ensure_dir(relnm):
    """ Accept relative filepath string, create it if it doesnt already exist
        return filepath string

    Args:
        relnm (str) : Relative name/path

    Returns:
        relnm (str)

    """

    d = os.path.join(os.getcwd(), relnm)
    if not os.path.exists(d):
        os.makedirs(d)
    return relnm


def check_list(values):
    values_str = '_'.join(values) if type(values) is list else values
    return values_str


def run_causal_path(relnm):
    print('Running CausalPath')
    out_f = open(os.path.join(relnm,'causal-path.out'),'w')
    args = ['java','-jar','/Users/estabroj/Aux_GIT_repos/causalpath/target/causalpath.jar', relnm]
    subprocess.run(args=args, stdout=out_f)
    out_f.close()
    print('Complete\n\n')


def loop_run_causal():
    contrasts = [map(str,list(t)) for t in zip(np.repeat(0, 15).tolist(),range(3, 48, 3))]
    contrasts_for_causal = ['significant-change-of-mean_{}_{}'.format(list(t)[0], list(t)[1]) for t in zip(np.repeat(0, 15).tolist(), range(3, 48, 3))]
    for contrast in contrasts:
        print('Running run_causal_series.py')
        args = ['python','run_causal_series.py','-vt','significant-change-of-mean','-tx',contrast[0],contrast[1]]
        subprocess.call(args=args)
        print('Complete\n\n')


def plot_pca():
    import seaborn as sns
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from sklearn.impute import SimpleImputer
    imp_med = SimpleImputer(missing_values = np.nan, strategy = 'median')

    pca = PCA(n_components = 2)

    meta = pd.read_csv('phospho_beta_meta.txt', sep = '\t', index_col = 0)
    phospho_data = load_phospho_data_raw()
    # phospho_data = load_phospho_data()

    phospho_filt = phospho_data.loc[:, meta.index]

    X = pca.fit_transform(imp_med.fit_transform(phospho_filt.T))
    pca_frame = pd.DataFrame(X, index = phospho_filt.T.index, columns = ['PCA-1', 'PCA-2'])
    joined_pca_frame = pca_frame.join(meta)

    plt.clf()

    ax = sns.scatterplot('PCA-1', 'PCA-2', data = joined_pca_frame, hue = 'Condition', style = 'Plex')
    plt.xlabel('{:10.4f}%'.format(pca.explained_variance_ratio_[0] * 100.0),size=5)
    plt.ylabel('{:10.4f}%'.format(pca.explained_variance_ratio_[1] * 100.0),size=5)
    plt.tight_layout()
    plt.legend(numpoints = 1,bbox_to_anchor=(1.1,1.), bbox_transform=plt.gcf().transFigure,fancybox = True, shadow = False,
               title = 'Categories', prop = {'size': 6}, title_fontsize = 4,
               frameon = False)

    plt.savefig('all_phospho_total_pca.png', dpi = 500, bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0., bbox_inches = 'tight')
    plt.clf()


def plot_pca_prot():
    import seaborn as sns
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from sklearn.impute import SimpleImputer
    imp_med = SimpleImputer(missing_values = np.nan, strategy = 'median')

    pca = PCA(n_components = 2)

    meta = pd.read_csv('phospho_beta_meta.txt', sep = '\t', index_col = 0)
    # phospho_data = load_total_protein()
    phospho_data = load_total_protein_raw()
    prot_map = {'WT_C1908_rep1':'plex1_WT_C1908','WT_C1908_rep2':'plex2_WT_C1908'}
    meta = meta.rename(prot_map)
    phospho_filt = phospho_data.loc[:, meta.index]

    X = pca.fit_transform(imp_med.fit_transform(phospho_filt.T))
    pca_frame = pd.DataFrame(X, index = phospho_filt.T.index, columns = ['PCA-1', 'PCA-2'])

    joined_pca_frame = pca_frame.join(meta)

    plt.clf()

    ax = sns.scatterplot('PCA-1', 'PCA-2', data = joined_pca_frame, hue = 'Type', style = 'Time')
    plt.xlabel('{:10.4f}%'.format(pca.explained_variance_ratio_[0] * 100.0),size=5)
    plt.ylabel('{:10.4f}%'.format(pca.explained_variance_ratio_[1] * 100.0),size=5)
    plt.tight_layout()
    plt.legend(numpoints = 1,bbox_to_anchor=(1.1,1.), bbox_transform=plt.gcf().transFigure,fancybox = True, shadow = False,
               title = 'Categories', prop = {'size': 6}, title_fontsize = 4,
               frameon = False)

    plt.savefig('all_prot_total_pca.png', dpi = 500, bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0., bbox_inches = 'tight')
    plt.clf()



def plot_pca_total_prot_phos():
    import seaborn as sns
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from sklearn.impute import SimpleImputer
    imp_med = SimpleImputer(missing_values = np.nan, strategy = 'median')

    pca = PCA(n_components = 2)

    meta = pd.read_csv('phospho_beta_meta.txt', sep = '\t', index_col = 0)

    # phospho_data = load_phospho_data_raw()
    phospho_data = pd.read_csv('phospho_prot_robust_z_new.txt', sep = '\t',
                               index_col = ['ID', 'Symbols', 'Sites', 'Effect'])
    phospho_filt = phospho_data.loc[:, meta.index]


    # prot_data = load_total_protein_raw()
    # prot_map = {'WT_C1908_rep1':'plex1_WT_C1908','WT_C1908_rep2':'plex2_WT_C1908'}
    # meta = meta.rename(prot_map)
    # prot_filt = prot_data.loc[:, meta.index]
    # phos_filt = phospho_filt.T.rename(prot_map).T
    #
    # phospho_filt = pd.concat([prot_filt,phos_filt])

    X = pca.fit_transform(imp_med.fit_transform(phospho_filt.T))
    pca_frame = pd.DataFrame(X, index = phospho_filt.T.index, columns = ['PCA-1', 'PCA-2'])

    joined_pca_frame = pca_frame.join(meta)

    plt.clf()

    ax = sns.scatterplot('PCA-1', 'PCA-2', data = joined_pca_frame, hue = 'Type', style = 'Time')
    plt.xlabel('{:10.4f}%'.format(pca.explained_variance_ratio_[0] * 100.0),size=5)
    plt.ylabel('{:10.4f}%'.format(pca.explained_variance_ratio_[1] * 100.0),size=5)
    plt.tight_layout()
    plt.legend(numpoints = 1,bbox_to_anchor=(1.1,1.), bbox_transform=plt.gcf().transFigure,fancybox = True, shadow = False,
               title = 'Categories', prop = {'size': 6}, title_fontsize = 4,
               frameon = False)

    plt.savefig('all_prot_phos_total_pca.png', dpi = 500, bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0., bbox_inches = 'tight')
    plt.clf()


def plot_pca_ortho_prot_phos():
    import seaborn as sns
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from sklearn.impute import SimpleImputer
    imp_med = SimpleImputer(missing_values = np.nan, strategy = 'median')

    pca = PCA(n_components = 2)

    meta = pd.read_csv('phospho_beta_meta.txt', sep = '\t', index_col = 0)

    phospho_filt = pd.read_csv('phospho_prot_data.txt',sep='\t')
    phospho_filt = phospho_filt.loc[:,meta.index]

    X = pca.fit_transform(imp_med.fit_transform(phospho_filt.T))
    pca_frame = pd.DataFrame(X, index = phospho_filt.T.index, columns = ['PCA-1', 'PCA-2'])

    joined_pca_frame = pca_frame.join(meta)

    plt.clf()

    ax = sns.scatterplot('PCA-1', 'PCA-2', data = joined_pca_frame, hue = 'Type', style = 'Time')
    plt.xlabel('{:10.4f}%'.format(pca.explained_variance_ratio_[0] * 100.0),size=5)
    plt.ylabel('{:10.4f}%'.format(pca.explained_variance_ratio_[1] * 100.0),size=5)
    plt.tight_layout()
    plt.legend(numpoints = 1,bbox_to_anchor=(1.1,1.), bbox_transform=plt.gcf().transFigure,fancybox = True, shadow = False,
               title = 'Categories', prop = {'size': 6}, title_fontsize = 4,
               frameon = False)

    plt.savefig('all_prot_phos_ortho_robust_z_pca.png', dpi = 500, bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0., bbox_inches = 'tight')
    plt.clf()



def generate_pallete_map(meta):
    groups = ['Time','Type']
    g_set = set()
    [g_set.update(meta[x].unique().tolist()) for x in groups]
    g_list = list(g_set)
    pal = sns.color_palette('colorblind',len(g_list))
    lut = dict(zip(map(str, g_list),pal))
    row_colors = []
    for g in groups:
        row_colors.append(meta[g].map(lut))

    row_frame = pd.concat(row_colors, axis = 1)

    return row_frame, lut


def generate_heatmap(phospho_filt, meta):
    out = 'total_phospho_heatmap.png'
    meta_filt = meta.reindex(phospho_filt.columns)

    imp_med = SimpleImputer(missing_values = np.nan, strategy = 'median')

    joined_frame, lut = generate_pallete_map(meta_filt)
    sns.set()

    filt_enrichment = pd.DataFrame(imp_med.fit_transform(phospho_filt),index=phospho_filt.index,columns = phospho_filt.columns).T
    ax = sns.clustermap(filt_enrichment, row_colors=joined_frame, cmap='coolwarm',cbar_kws={'label': 'Phosphosites'},**{'yticklabels':False,'xticklabels':False,'center':0.0}) #, method=method, metric=metric)
    font = mpl.font_manager.FontProperties(family = 'times new roman', style = 'italic', size = 5)
    text = ax.fig.axes[-1].yaxis.label
    text.set_font_properties(font)
    plt.setp(ax.ax_row_colors.get_xticklabels(), size = 5)
    plt.setp(ax.ax_heatmap.get_xticklabels(), rotation = 90, size=5)
    plt.setp(ax.ax_heatmap.get_yticklabels(), size=5)
    markers = [plt.Line2D([0, 0], [0, 0], color = color, marker = 's', linestyle = '') for color in lut.values()]
    plt.legend(markers, lut.keys(), numpoints = 1,bbox_to_anchor=(1.,1.), bbox_transform=plt.gcf().transFigure,fancybox = False, shadow = False,
               title = 'Categories', ncol = 4, prop = {'size': 6}, title_fontsize = 4,
               frameon = False)

    fig = mpl.pyplot.gcf()
    fig.set_size_inches(14.5, 10.5)
    ax.savefig(out,dpi=500)



def generate_heatmap_prot(meta):
    prot_data = load_total_protein_raw()
    prot_map = {'WT_C1908_rep1':'plex1_WT_C1908','WT_C1908_rep2':'plex2_WT_C1908'}
    meta = meta.rename(prot_map)
    phospho_filt = prot_data.loc[:, meta.index]

    out = 'total_prot_heatmap.png'
    meta_filt = meta.reindex(phospho_filt.columns)

    imp_med = SimpleImputer(missing_values = np.nan, strategy = 'median')

    joined_frame, lut = generate_pallete_map(meta_filt)
    sns.set()

    filt_enrichment = pd.DataFrame(imp_med.fit_transform(phospho_filt),index=phospho_filt.index,columns = phospho_filt.columns).T
    ax = sns.clustermap(filt_enrichment, z_score=1, row_colors=joined_frame, cmap='coolwarm',cbar_kws={'label': 'Protein'},**{'yticklabels':False,'xticklabels':False}) #, method=method, metric=metric)
    font = mpl.font_manager.FontProperties(family = 'times new roman', style = 'italic', size = 5)
    text = ax.fig.axes[-1].yaxis.label
    text.set_font_properties(font)
    plt.setp(ax.ax_row_colors.get_xticklabels(), size = 5)
    plt.setp(ax.ax_heatmap.get_xticklabels(), rotation = 90, size=5)
    plt.setp(ax.ax_heatmap.get_yticklabels(), size=5)
    markers = [plt.Line2D([0, 0], [0, 0], color = color, marker = 's', linestyle = '') for color in lut.values()]
    plt.legend(markers, lut.keys(), numpoints = 1,bbox_to_anchor=(1.,1.), bbox_transform=plt.gcf().transFigure,fancybox = False, shadow = False,
               title = 'Categories', ncol = 4, prop = {'size': 6}, title_fontsize = 4,
               frameon = False)

    fig = mpl.pyplot.gcf()
    fig.set_size_inches(14.5, 10.5)
    ax.savefig(out,dpi=500)



def plot_peptide(target):
    target_frame = data_stack_joined.loc[target]
    sns.catplot(x = 'Type', y = 'Value', col = 'Time', data = target_frame)
    plt.show()

def corrfunc(x, y, **kws):
    r, _ = stats.pearsonr(x, y)
    ax = plt.gca()
    ax.annotate("r = {:.2f}".format(r),
                xy=(.1, .9), xycoords=ax.transAxes)

def map_color(meta):
    meta_dict = meta.Condition.to_dict()
    pal = sns.color_palette('colorblind',len(meta_dict.values()))
    color_map = dict(zip(meta_dict.keys(),pal))

    return color_map


def plot_pairplot(df, out, meta):
    replacements = meta.Condition.to_dict()

    plt.clf()
    g = sns.PairGrid(df, palette = ["red"])
    g.map_upper(plt.scatter, s = 5, alpha=0.5, edgecolor='k')
    g.map_diag(sns.distplot, kde = True)
    g.map_lower(sns.kdeplot, cmap = "Blues_d", n_levels=50, **{"linewidths":0.1})
    g.map_lower(corrfunc)

    for i in range(len(replacements)):
        for j in range(len(replacements)):
            xlabel = g.axes[i][j].get_xlabel()
            ylabel = g.axes[i][j].get_ylabel()
            if xlabel in replacements.keys():
                g.axes[i][j].set_xlabel(replacements[xlabel])
            if ylabel in replacements.keys():
                g.axes[i][j].set_ylabel(replacements[ylabel])


    g.savefig(out, dpi = 500)


def draw_directed_labeled_edge_graph(shadow_pair_dict, total_enrichment):
    increase_edge, decrease_edge, synergy_edge, null_edge = bin_edge_types(shadow_pair_dict)
    graph = nx.DiGraph()

    graph.add_edges_from(increase_edge, label = 'I')
    graph.add_edges_from(decrease_edge, label = 'I')
    graph.add_edges_from(synergy_edge, label = 'S')
    # graph.add_edges_from(null_edge, label='N')
    pos = nx.nx_agraph.graphviz_layout(graph, prog = 'neato')

    labels = {node: node for node in graph.nodes()}
    nx.set_node_attributes(graph, 'labels', labels.values())
    nx.draw_networkx_labels(graph, pos, labels, font_size = 8, font_color='grey',font_weight='bold')
    val_map = total_enrichment.loc[:, list(labels.keys())].median().to_dict()
    values = [val_map.get(node, 0.0) for node in graph.nodes()]
    edge_labels = dict([((u, v,), d['label']) for u, v, d in graph.edges(data = True)])
    edge_list = [(u, v) for u, v in graph.edges()]

    nx.draw_networkx_edge_labels(graph, pos, edge_labels = edge_labels, font_size = 3)

    nodes = nx.draw_networkx_nodes(graph, pos, cmap = plt.get_cmap('plasma'), node_color = values, node_size = 500,
                                   alpha = 1, with_labels = True)

    i_idx = ['_'.join(x) for x in increase_edge]
    d_idx = ['{}_{}'.format(x[1], x[0]) for x in decrease_edge]
    s_idx = ['{}_{}'.format(x[1], x[0]) for x in synergy_edge]
    # s_idx = ['_'.join(x) for x in synergy_edge]

    idxs_r = i_idx + d_idx + s_idx

    colors = filt_shadow.loc[idxs_r]['correlation'].values.tolist()

    i_colors = filt_shadow.loc[i_idx]['correlation'].values.tolist()
    d_colors = filt_shadow.loc[d_idx]['correlation'].values.tolist()
    s_colors = filt_shadow.loc[s_idx]['correlation'].values.tolist()


    edges = nx.drawing.nx_pylab.draw_networkx_edges(graph, pos, node_size = 300, edge_color = colors,
                                                    edgelist = edge_list, width = 0.01, edge_cmap = plt.cm.coolwarm,
                                                    arrowstyle = '-', arrows = False, edge_vmin = min(colors),
                                                    edge_vmax = max(colors))

    # increase_collection = nx.drawing.nx_pylab.draw_networkx_edges(graph, pos, node_size = 300, edgelist = increase_edge,
    #                                                               edge_color = i_colors, width = 2,
    #                                                               edge_cmap = plt.cm.coolwarm, edge_vmin = min(colors),
    #                                                               edge_vmax = max(colors), arrowstyle = '-|>',
    #                                                               style = 'dotted')
    # for patch in increase_collection:
    #     patch.set_linestyle('dotted')

    decrease_collection = nx.drawing.nx_pylab.draw_networkx_edges(graph, pos, node_size = 300, edgelist = decrease_edge,
                                                                  edge_color = d_colors, width = 2,
                                                                  edge_cmap = plt.cm.coolwarm, edge_vmin = min(colors),
                                                                  edge_vmax = max(colors), arrowstyle = '-|>',
                                                                  style = 'dashed')
    for patch in decrease_collection:
        patch.set_linestyle('dashed')

    nx.drawing.nx_pylab.draw_networkx_edges(graph, pos, node_size = 300, edgelist = synergy_edge, edge_color = s_colors,
                                            width = 2, edge_cmap = plt.cm.coolwarm, edge_vmin = min(colors),
                                            edge_vmax = max(colors), arrowstyle = '<|-|>')


    edge_cb = plt.colorbar(edges)
    edge_cb.set_label('Weight correlation')

    node_cb = plt.colorbar(nodes)
    node_cb.set_label('Node median RES')
    plt.axis('off')



def build_example_shadow():


    regulator_a = [('Regulator-A', 'gene 1'), ('Regulator-A', 'gene 2'), ('Regulator-A', 'gene 3'),
                   ('gene 3', 'gene 4'),('gene 3','gene 5')]

    nodelist_i = ['gene 1', 'gene 2','gene 3']
    nodelist_u = ['gene 4', 'gene 5']
    main_nodes = ['Regulator-A']

    activated = [ ('Regulator-A', 'gene 3'),
                   ('gene 3', 'gene 4'),('gene 3','gene 5')]
    repressed = [('Regulator-A', 'gene 1')]
    non_mon = [('Regulator-A', 'gene 2')]

    graph = nx.DiGraph()
    graph.add_edges_from(activated, label = 'A', color='r')
    graph.add_edges_from(repressed, label = 'R', color='blue')
    graph.add_edges_from(non_mon, label = '-',color='grey')
    pos = nx.nx_agraph.graphviz_layout(graph, prog = 'neato')

    labels = {node: node for node in graph.nodes()}
    nx.set_node_attributes(graph, 'labels', labels.values())
    nx.draw_networkx_labels(graph, pos, labels, font_size = 7)

    nodes = nx.draw_networkx_nodes(graph, pos, nodelist = nodelist_i, node_size = 2000, alpha = 1, with_labels = True,
                                   node_color = 'grey')
    nodes = nx.draw_networkx_nodes(graph, pos, nodelist = nodelist_u, node_size = 2000, alpha = 1, with_labels = True,
                                   node_color = 'grey')
    nodes = nx.draw_networkx_nodes(graph, pos, nodelist = main_nodes, node_size = 2000, alpha = 1, with_labels = True,
                                   node_color = 'red')

    increase_collection = nx.draw_networkx_edges(graph, pos, node_size = 2000, edgelist = activated,
                                                                  edge_color = ['red','red','red'], width = 2, arrowstyle = '-|>')
    decrease_collection = nx.draw_networkx_edges(graph, pos, node_size = 2000, edgelist = repressed,
                                                                  edge_color = ['blue'], width = 2, arrowstyle = '-[')
    increase_collection = nx.draw_networkx_edges(graph, pos, node_size = 2000, edgelist = non_mon,
                                                                  edge_color = 'grey', width = 2, arrowstyle = '-')

    plt.axis('off')


def plot_moa_kde():
    nsamples = 10000
    means = [-.8]
    sds = [.22]
    weights = [0.05]
    draws = np.random.multinomial(nsamples, weights)
    moa_down = np.concatenate(list(starmap(np.random.normal, zip(means, sds, draws))))

    means = [0]
    sds = [.2]
    weights = [0.8]
    draws = np.random.multinomial(nsamples, weights)
    moa_non = np.concatenate(list(starmap(np.random.normal, zip(means, sds, draws))))

    means = [.9]
    sds = [.25]
    weights = [0.3]
    draws = np.random.multinomial(nsamples, weights)
    moa_up = np.concatenate(list(starmap(np.random.normal, zip(means, sds, draws))))

    fig, ax = plt.subplots(figsize = (7, 14))

    ax = sns.kdeplot(moa_down, color = 'b', shade=True,label='Repressed')
    ax = sns.kdeplot(moa_non, color = 'grey', shade=True,label='Non-monotonically regulated')
    ax = sns.kdeplot(moa_up, color = 'r', shade=True, label='Activated')
    sns.despine()
    fig.savefig('example_moa.pdf', dpi = 300,bbox_inches = 'tight')



def build_example_shadow():
    graph = nx.DiGraph()

    regulator_a = [('Regulator-A', 'gene 1'), ('Regulator-A', 'gene 2'), ('Regulator-A', 'gene 3'),
                   ('Regulator-A', 'gene 4')]
    regulator_b = [('Regulator-B', 'gene 1'), ('Regulator-B', 'gene 3'), ('Regulator-B', 'gene 5'),
                   ('Regulator-B', 'gene 6')]

    nodelist_i = ['gene 1', 'gene 3']
    nodelist_u = ['gene 4', 'gene 2', 'gene 5', 'gene 6']
    main_nodes = ['Regulator-B', 'Regulator-A']

    graph.add_edges_from(regulator_a, label = 'A')
    graph.add_edges_from(regulator_b, label = 'B')
    pos = nx.nx_agraph.graphviz_layout(graph, prog = 'neato')
    labels = {node: node for node in graph.nodes()}
    nx.set_node_attributes(graph, 'labels', labels.values())
    nx.draw_networkx_labels(graph, pos, labels, font_size = 7)

    edge_labels = dict([((u, v,), d['label']) for u, v, d in graph.edges(data = True)])
    edge_list = [(u, v) for u, v in graph.edges()]

    nx.draw_networkx_edge_labels(graph, pos, edge_labels = edge_labels, font_size = 3)

    nodes = nx.draw_networkx_nodes(graph, pos, nodelist = nodelist_i, node_size = 2000, alpha = 1, with_labels = True,
                                   node_color = 'green')
    nodes = nx.draw_networkx_nodes(graph, pos, nodelist = nodelist_u, node_size = 2000, alpha = 1, with_labels = True,
                                   node_color = 'grey')
    nodes = nx.draw_networkx_nodes(graph, pos, nodelist = main_nodes, node_size = 2000, alpha = 1, with_labels = True,
                                   node_color = 'red')

    increase_collection = nx.drawing.nx_pylab.draw_networkx_edges(graph, pos, node_size = 2000, edgelist = regulator_a,
                                                                  edge_color = 'black', width = 2, arrowstyle = '-|>')

    decrease_collection = nx.drawing.nx_pylab.draw_networkx_edges(graph, pos, node_size = 2000, edgelist = regulator_b,
                                                                  width = 2, edge_color = 'black', arrowstyle = '-|>')

    plt.axis('off')
    targets = ['MAML1', 'BATF3', 'GATA6', 'ARNT', 'EPAS1', 'MED12', 'MED21', 'MAP2K4']

    regul = read_pickle('DCIS_bolstered_regulon.pkl')
    regul = regul.reset_index()
    total_graph = nx.from_pandas_edgelist(regul,'Regulator','Target',edge_attr='MoA',create_using=nx.DiGraph())
    only_containing_nodes = lambda x: 'MAML1' in x and 'BATF3' in x and 'GATA6' and 'ARNT' in x and 'EPAS1' in x and 'MED12' in x and 'MED21' in x and 'MAP2K4' in x
    G = total_graph
    all_simple_paths = nx.all_simple_paths(G, source='EPAS1',target='MAML1')
    all_shortest_paths = nx.all_shortest_paths(G, source='EPAS1',target='MAML1')


def build_example_shadow_repressive():
    graph = nx.DiGraph()

    regulator_a = [('Regulator-A', 'gene 2'), ('Regulator-A', 'gene 4')]
    regulator_b = [('Regulator-B', 'gene 5'), ('Regulator-B', 'gene 6')]

    nodelist_a = ['Regulator-B']
    nodelist_b = ['Regulator-A']
    nodelist_u = ['gene 4', 'gene 2', 'gene 5', 'gene 6']

    graph.add_edges_from(regulator_a, label = 'A')
    graph.add_edges_from(regulator_b, label = 'B')
    pos = nx.nx_agraph.graphviz_layout(graph, prog = 'neato')
    labels = {node: node for node in graph.nodes()}
    labels.pop('gene 1')
    labels.pop('gene 3')
    nx.set_node_attributes(graph, 'labels', labels.values())
    nx.draw_networkx_labels(graph, pos, labels, font_size = 7)

    edge_labels = dict([((u, v,), d['label']) for u, v, d in graph.edges(data = True)])
    edge_labels.pop(('Regulator-A', 'gene 1'))
    edge_labels.pop(('Regulator-A', 'gene 3'))
    edge_labels.pop(('Regulator-B', 'gene 3'))
    edge_labels.pop(('Regulator-B', 'gene 1'))

    nx.draw_networkx_edge_labels(graph, pos, edge_labels = edge_labels, font_size = 3)

    nodes = nx.draw_networkx_nodes(graph, pos, nodelist = nodelist_a, node_size = 2000, alpha = 1, with_labels = True,
                                   node_color = 'r')

    nodes = nx.draw_networkx_nodes(graph, pos, nodelist = nodelist_b, node_size = 2000, alpha = 1, with_labels = True,
                                   node_color = 'grey')

    nodes = nx.draw_networkx_nodes(graph, pos, nodelist = nodelist_u, node_size = 2000, alpha = 1, with_labels = True,
                                   node_color = 'grey')

    increase_collection = nx.drawing.nx_pylab.draw_networkx_edges(graph, pos, node_size = 2000, edgelist = regulator_a,
                                                                  edge_color = 'black', width = 2, arrowstyle = '-|>')

    decrease_collection = nx.drawing.nx_pylab.draw_networkx_edges(graph, pos, node_size = 2000, edgelist = regulator_b,
                                                                  width = 2, edge_color = 'black', arrowstyle = '-|>')

    plt.axis('off')


def plot_causal(results_f, regulator, value_changes):
    import pandas as pd
    import networkx as nx
    import matplotlib.pyplot as plt

    results = pd.read_csv(results_f, sep = '\t')
    total_graph = nx.from_pandas_edgelist(results, 'Source', 'Target',
                                          edge_attr = ['Relation', 'Sites', 'Source data ID', 'Source change',
                                                       ' Source change pval', 'Target data ID', 'Target change',
                                                       'Target change pval'], create_using = nx.DiGraph())
    val_map = value_changes['Change amount'].to_dict()
    graph = total_graph.subgraph(nx.shortest_path(total_graph.to_undirected(), regulator))
    graph_copy = graph.copy()

    # pop_nodes = [x for x in graph_degree if graph_degree[x] <= 4]
    # [graph_copy.remove_node(x) for x in pop_nodes]
    #  pos = nx.nx_agraph.graphviz_layout(graph_copy, prog = 'sfdp', args = '-Goverlap=false')
    # pos = nx.nx_agraph.graphviz_layout(graph_copy)#, prog = 'sfdp', args = '-Goverlap=false')
    pos = nx.spring_layout(graph_copy, scale=1000)#, prog = 'sfdp', args = '-Goverlap=false')
    # pos = nx.nx_agraph.graphviz_layout(graph_copy, args = '-Goverlap=false len=2.0') #, prog = 'sfdp', args = '-Goverlap=false')
    values = [val_map.get(node, 0.0) for node in graph_copy.nodes()]
    labels = {node: node for node in graph_copy.nodes()}

    nx.set_node_attributes(graph_copy, 'labels', labels.values())
    # nodes = nx.draw_networkx_nodes(graph_copy, pos, node_size = 400, alpha = 1, with_labels = True)
    nodes = nx.draw_networkx_nodes(graph, pos, cmap = plt.get_cmap('plasma'), node_color = values, node_size = 300,
                                   alpha = 1, with_labels = True)



    nx.draw_networkx_labels(graph_copy, pos, labels, font_size = 7)
    colors = list(nx.get_edge_attributes(graph_copy, 'Source change').values())
    edges = nx.drawing.nx_pylab.draw_networkx_edges(graph_copy, pos, edge_color = colors, width = 2,
                                                    edge_cmap = plt.cm.coolwarm, arrowstyle = '-|>', arrows = False)

    edge_cb = plt.colorbar(edges)
    edge_cb.set_label('Weight correlation')

    node_cb = plt.colorbar(nodes)
    node_cb.set_label('Node median RES')

    plt.axis('off')


def split_fit_predictions(expr, total_enrichment, regulon, regulator, n_splits = 5, n_repeats = 10, regressor = 'bayesridge'):
    """

    Args:
        expr (:obj: `pandas DataFrame`): pandas DataFrame of shape [n_samps, n_feats]
        zframe (:obj: `pandas DataFrame`): pandas DataFrame of activity scores for specified regulator
        regulon (:obj: `pandas DataFrame`): pandas DataFrame of regulon returned by compile_regulon
            with columns ['Target', 'MoA', 'likelihood']
        n_splits (int) : Number of splits for each cross-fold validation
        n_repeats (int) : Number of repeated cross-fold validations
        regressor (str) : sklearn regressor used to fit and predict NES

    Returns:
        scores (np.array) : np.array of length [n_splits * n_repeats] using default RegressorMixin score (R2)
        deviation (float) : standard deviation of scores

    """
    zframe = total_enrichment.loc[:,regulator]
    sub_reg = regulon.loc[regulator]

    regressor_opt = {'ridge': linear_model.Ridge(), 'bayesridge': linear_model.BayesianRidge(),
                     'ols': linear_model.LinearRegression(), 'lasso': linear_model.Lasso(),
                     'RandomForest': RandomForestRegressor(max_depth=2, random_state=0, n_estimators=100),
                     'SGD':SGDRegressor(max_iter=50, penalty=None, eta0=0.1)}
    l_rg = regressor_opt[regressor]

    y = zframe
    X = expr.loc[:, sub_reg.Target].dropna(axis = 1)
    l_rg.fit(X, y)
    rkf = RepeatedKFold(n_splits = n_splits, n_repeats = n_repeats)
    kf = KFold(n_splits, shuffle = True)

    predictions = []
    for i in range(n_repeats):
        pred = pd.DataFrame(cross_val_predict(l_rg, X, y, cv = kf),index=y.index,columns = ['{}'.format(i)])
        predictions.append(pred)

    pred_frame = pd.concat(predictions, axis=1)
    scores = cross_val_score(l_rg, X, y, cv = rkf)

    return scores, pred_frame



from matplotlib import cbook
from matplotlib.colors import Normalize
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.offsetbox
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import AxesGrid


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap



class MidPointNorm(Normalize):
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self,vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = np.getmask(result)
                result = np.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            #First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint
            resdat[resdat>0] /= abs(vmax - midpoint)
            resdat[resdat<0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = np.array(resdat, mask=result.mask, copy=False)

        if is_scalar:
            result = result[0]
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if cbook.iterable(value):
            val = np.asarray(value)
            val = 2 * (val-0.5)
            val[val>0]  *= abs(vmax - midpoint)
            val[val<0] *= abs(vmin - midpoint)
            val += midpoint
            return val
        else:
            val = 2 * (val - 0.5)
            if val < 0:
                return  val*abs(vmin-midpoint) + midpoint
            else:
                return  val*abs(vmax-midpoint) + midpoint


class AnchoredHScaleBar(matplotlib.offsetbox.AnchoredOffsetbox):
    """ size: length of bar in data units
        extent : height of bar ends in axes units """
    def __init__(self, size=1, extent = 0.03, label="", loc=2, ax=None,
                 pad=0.4, borderpad=0.5, ppad = 0, sep=2, prop=None,
                 frameon=True, **kwargs):
        if not ax:
            ax = plt.gca()
        trans = ax.get_xaxis_transform()
        size_bar = matplotlib.offsetbox.AuxTransformBox(trans)
        line = Line2D([0,size],[0,0], **kwargs)
        vline1 = Line2D([0,0],[-extent/2.,extent/2.], **kwargs)
        vline2 = Line2D([size,size],[-extent/2.,extent/2.], **kwargs)
        size_bar.add_artist(line)
        size_bar.add_artist(vline1)
        size_bar.add_artist(vline2)
        txt = matplotlib.offsetbox.TextArea(label, minimumdescent=False)
        self.vpac = matplotlib.offsetbox.VPacker(children=[size_bar,txt],
                                 align="center", pad=ppad, sep=sep)
        matplotlib.offsetbox.AnchoredOffsetbox.__init__(self, loc, pad=pad,
                 borderpad=borderpad, child=self.vpac, prop=prop, frameon=frameon)


def compare_correlations(total_enrichment, log_expr, r, meta):
    import matplotlib.patches as mpatches

    infer_vals = total_enrichment[r]
    meta_filt = meta[meta.RNA_SeqID.isin(infer_vals.index)]
    meta_filt.set_index('RNA_SeqID', inplace = True)
    meta_filt = meta_filt.reindex(infer_vals.index)
    infer_vals = pd.concat([infer_vals, log_expr.loc[r, :]], axis = 1)
    infer_vals.columns = ['{}-nes'.format(r), '{}-expr'.format(r)]

    joined_infer_frame = pd.concat([infer_vals, meta_filt['condition']], axis = 1)

    translate_dict = {'ctd2': 'patient', 'cd34_pooled_technical': 'CD34+p', 'Individual_MNC': 'MNC', 'Individual_CD34': 'CD34'}
    color_dict = {'ctd2': "Greys", 'cd34_pooled_technical': "Blues", 'Individual_MNC': "Greens",
                      'Individual_CD34': 'red'}
    fig, axs = plt.subplots(2,1, figsize=(5,8))
    label_patches = []
    for group,frame in joined_infer_frame.groupby('condition'):
        trans_frame = frame.replace({'condition':translate_dict})
        c = color_dict[group]
        if group != 'Individual_CD34':
            ax = sns.kdeplot(trans_frame['{}-nes'.format(r)], trans_frame['{}-expr'.format(r)], cmap = color_dict[group]+'_d', shade = False, shade_lowest = False,ax = axs[0],**{"linewidths":0.5})
            label_patch = mpatches.Patch(color = sns.color_palette(color_dict[group])[2], label = translate_dict[group])
            label_patches.append(label_patch)
        else:
            ax = sns.scatterplot(x = '{}-nes'.format(r), y = '{}-expr'.format(r), data = trans_frame, color = color_dict[group], ax=axs[0], size=5, markers ='*')
            label_patch = mpatches.Patch(color = color_dict[group], label = translate_dict[group])
            label_patches.append(label_patch)
    ax.legend(handles = label_patches, loc = 'center left', bbox_to_anchor=(1,0.5),prop={'size':5})
    ax.set_xlabel('RES',fontsize=6)
    ax.set_ylabel('log2 TPM',fontsize=6)
    plt.setp(axs[0].xaxis.get_majorticklabels(),size=6)
    plt.setp(axs[0].yaxis.get_majorticklabels(),size=6)

    ax2 = sns.swarmplot(y = '{}-nes'.format(r), x = 'condition'.format(r), data=joined_infer_frame.replace({'condition':translate_dict}), size=3, ax=axs[1],palette={'patient': 'grey', 'CD34+p' : 'skyblue', 'MNC' : 'seagreen','CD34':'red'})
    infer_vals = total_enrichment[r]
    meta_filt = meta[meta.RNA_SeqID.isin(infer_vals.index)]
    meta_filt.set_index('RNA_SeqID', inplace = True)
    meta_filt = meta_filt.reindex(infer_vals.index)
    infer_vals = pd.concat([infer_vals, log_expr.loc[r, :]], axis = 1)
    infer_vals.columns = ['{}-nes'.format(r), '{}-expr'.format(r)]

    joined_infer_frame = pd.concat([infer_vals, meta_filt['condition']], axis = 1)

    translate_dict = {'ctd2': 'patient', 'cd34_pooled_technical': 'CD34+p', 'Individual_MNC': 'MNC', 'Individual_CD34': 'CD34'}
    color_dict = {'ctd2': "Greys", 'cd34_pooled_technical': "Blues", 'Individual_MNC': "Greens",
                      'Individual_CD34': 'red'}
    fig, axs = plt.subplots(2,1, figsize=(5,8))
    label_patches = []
    for group,frame in joined_infer_frame.groupby('condition'):
        trans_frame = frame.replace({'condition':translate_dict})
        c = color_dict[group]
        if group != 'Individual_CD34':
            ax = sns.kdeplot(trans_frame['{}-nes'.format(r)], trans_frame['{}-expr'.format(r)], cmap = color_dict[group]+'_d', shade = False, shade_lowest = False,ax = axs[0],**{"linewidths":0.5})
            label_patch = mpatches.Patch(color = sns.color_palette(color_dict[group])[2], label = translate_dict[group])
            label_patches.append(label_patch)
        else:
            ax = sns.scatterplot(x = '{}-nes'.format(r), y = '{}-expr'.format(r), data = trans_frame, color = color_dict[group], ax=axs[0], size=5, markers ='*')
            label_patch = mpatches.Patch(color = color_dict[group], label = translate_dict[group])
            label_patches.append(label_patch)
    ax.legend(handles = label_patches, loc = 'center left', bbox_to_anchor=(1,0.5),prop={'size':5})
    ax.set_xlabel('RES',fontsize=6)
    ax.set_ylabel('log2 TPM',fontsize=6)
    plt.setp(axs[0].xaxis.get_majorticklabels(),size=6)
    plt.setp(axs[0].yaxis.get_majorticklabels(),size=6)

    palette = {'patient': 'grey', 'CD34+p': 'skyblue', 'MNC': 'seagreen', 'CD34': 'red'}
    ax2 = sns.violinplot(y = '{}-nes'.format(r), x = 'condition', data = trans_frame, palette=palette, ax = axs[1])
    ax2 = sns.swarmplot(y = '{}-nes'.format(r), x = 'condition'.format(r), data=joined_infer_frame.replace({'condition':translate_dict}), size=3, ax=axs[1],palette={'patient': 'grey', 'CD34+p' : 'skyblue', 'MNC' : 'seagreen','CD34':'red'})
    ax2.set_xlabel('Condtiion',fontsize=6)
    ax2.set_ylabel('{} RES'.format(r),fontsize=6)

    plt.setp(axs[1].xaxis.get_majorticklabels(), rotation = 45, size=6)
    plt.setp(axs[1].yaxis.get_majorticklabels(),size=6)
    fig.tight_layout()
    ax2.set_xlabel('Condtiion',fontsize=6)
    ax2.set_ylabel('{} RES'.format(r),fontsize=6)

    plt.setp(axs[1].xaxis.get_majorticklabels(), rotation = 45, size=6)
    plt.setp(axs[1].yaxis.get_majorticklabels(),size=6)
    fig.tight_layout()

    fig.savefig('{}_swarm.pdf'.format(r))

    sns.violinplot(x = '{}-nes'.format(r), y = '{}-expr'.format(r), data = trans_frame,
                   palette = palette)


