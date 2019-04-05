from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
import matplotlib.pyplot as plt
from decimal import Decimal
from matplotlib.patches import Patch
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D

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
    # ax.text(ax.get_xlim()[1] * 0.91, plt_ymax * 0.82, size = 15, ha = 'right',
    #         s = "accuracy: {:2.3f} ".format(np.mean(scores)))

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
    path = '/Users/estabroj/PycharmProjects/regulon_enrichment/experiments/CTD2/data'
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
edges = nx.drawing.nx_pylab.draw_networkx_edges(graph, pos, edge_color = colors, width = 2, edge_cmap = plt.cm.coolwarm,
                                                arrowstyle = '-|>', arrows = False)

edge_cb = plt.colorbar(edges)
edge_cb.set_label('Weight correlation')

node_cb = plt.colorbar(nodes)
node_cb.set_label('Node median RES')
plt.axis('off')


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
    from sklearn import linear_model
    from sklearn.model_selection import KFold, cross_val_predict, RepeatedKFold, cross_val_score

    regressor_opt = {'ridge': linear_model.Ridge(), 'bayesridge': linear_model.BayesianRidge(),
                     'ols': linear_model.LinearRegression(), 'lasso': linear_model.Lasso()}
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