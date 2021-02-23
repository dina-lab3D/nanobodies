
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from statannot import add_stat_annotation


# no_cut = pd.read_csv("/cs/labs/dina/tomer.cohen13/nanobodies/COVID_19/nb_paper/epiDock_summey.csv", skipinitialspace=True)
# no_cut = no_cut.drop([str(i) for i in range(1,320)] + [str(i) for i in range(541,1001)] + ["Unnamed: 1002"], axis=1)
# no_cut.to_csv("/cs/labs/dina/tomer.cohen13/nanobodies/COVID_19/nb_paper/epiDock_summey_cut.csv")
# exit()

########################################################################################################################
#                                                   Consurf                                                            #
########################################################################################################################


consurf_df = pd.read_csv("/cs/labs/dina/tomer.cohen13/Epitop_paper/consurf.grades", skipfooter=5, skiprows=list(range(13)) + [14], delim_whitespace=True, skipinitialspace=True)

positions = []
consurf_scores = []

for score, pos in zip(consurf_df["SCORE"],consurf_df["3LATOM"]):
    if pos == '-':
        continue

    positions.append(int(pos[3:6]))
    consurf_scores.append(score)

########################################################################################################################
#                                                   Mutations                                                          #
########################################################################################################################

mutations_df = pd.read_csv("/cs/usr/tomer.cohen13/lab/Epitop_paper/coronavirus3d-S-mutations.tsv",sep='\t')
only_ms = mutations_df[mutations_df['Annotation Type'].map(lambda x: 'missense_variant' in x)]

mutations = np.zeros(len(positions), dtype=int)
start_pos = positions[0]
end_pos = positions[-1]

for row in only_ms.index:

    already_visited = set()
    row_mut = only_ms["AA change"][row].split(";")
    for mut in row_mut:
        pos = mut[13:].split("-")
        if len(pos) == 1:
            pos = int(''.join(c for c in pos[0] if c.isdigit()))
            if (start_pos <= pos <= end_pos) and (pos not in already_visited):
                mutations[pos-start_pos] += only_ms["Count"][row]
                already_visited.add(pos)


########################################################################################################################
#                                                  Ab Epitops                                                          #
########################################################################################################################


epitops_df = pd.read_csv("/cs/usr/tomer.cohen13/lab/Epitop_paper/epiDock_summery_cut.csv")
new_ab_nr = pd.read_csv("/cs/usr/tomer.cohen13/lab/Epitop_paper/new_epiDock_summey_cut.csv")
epitops_df['ANTIBODY'] = epitops_df['ANTIBODY'].map(lambda x: x.split('pdb')[1].split(".ent")[0].upper())
new_ab_nr['ANTIBODY'] = new_ab_nr['ANTIBODY'].map(lambda x: x.split('.')[0].upper())

epitops_df = pd.concat([epitops_df, new_ab_nr])

non_redundent = list(pd.read_csv("/cs/usr/tomer.cohen13/lab/Epitop_paper/antibodies_nr.csv")['PDB ID'])
epitops_nr = epitops_df[epitops_df['ANTIBODY'].map(lambda x: x in non_redundent)]
epitops_nr = epitops_nr.sort_values(by="ANTIBODY").reset_index().drop([22,10])  # drop 6zcz nanobody, and 7che 2

mean_epitop_nr = np.zeros(len(positions), dtype=float)
for i in range(len(positions)):
    mean_epitop_nr[i] = np.mean(epitops_nr[str(positions[i])])


########################################################################################################################
#                                                  Nb Epitops                                                          #
########################################################################################################################


nb_epitops_df = pd.read_csv("/cs/labs/dina/tomer.cohen13/Epitop_paper/nb_epiDock_summey_cut.csv").reindex()
nb_epitops_df = (nb_epitops_df.drop([12]))  # 7JWB
# epitops_df['ANTIBODY'] = epitops_df['ANTIBODY'].map(lambda x: x.split('pdb')[1].split(".ent")[0].upper())

mean_nb_epitop_nr = np.zeros(len(positions), dtype=float)
for i in range(len(positions)):
    mean_nb_epitop_nr[i] = np.mean(nb_epitops_df[str(positions[i])])



########################################################################################################################
#                                                   Ace2 binding                                                       #
########################################################################################################################

ace2_df_bind = pd.read_csv("/cs/labs/dina/tomer.cohen13/Epitop_paper/ace2_binding.csv")
ace2_bind = np.zeros(len(positions))
start_pos = positions[0]

for i in range(len(positions)):
    pos_df = ace2_df_bind[ace2_df_bind["site_SARS2"] == (i+start_pos)]
    max_bind = pos_df["bind_avg"].max()
    ace2_bind[i] = max_bind


########################################################################################################################
#                                                   Ace2 contact                                                       #
########################################################################################################################

ace2_contact_pos = [417,446,449,453,455,456,475,476,484,486,487,489,490,493,496,497,498,500,501,502,505]

ace2_contact = np.zeros(len(positions))
start_pos = positions[0]

for contact_pos in ace2_contact_pos:
    ace2_contact[contact_pos-start_pos] = 1


########################################################################################################################
#                                               curvature boxplot                                                      #
########################################################################################################################

ab_nr_data = pd.read_csv("/cs/usr/tomer.cohen13/lab/Epitop_paper/PDB_NR_new_ab.csv")
nb_nr_data = pd.read_csv("/cs/usr/tomer.cohen13/lab/Epitop_paper/PDB_NR_new_nb.csv")

nb_nr_data = (nb_nr_data.drop([7]))


box_data  = pd.DataFrame({"PDB": list(nb_nr_data["PDB"]) + list(ab_nr_data["PDB"]), "name": list(nb_nr_data["Nb name"]) + list(ab_nr_data["Ab name"]), "type": len(nb_nr_data) * ["Nb"] + len(ab_nr_data) * ["Fab"] , "Spike cavity": list(nb_nr_data["cavity spike"]) + list(ab_nr_data["cavity spike"]) })

ax = sns.boxplot(x="type", y="Spike cavity", data=box_data, color="white")
ax2 = sns.swarmplot(x="type", y="Spike cavity", data=box_data)
sns.despine(bottom = False, left = False, right=True, top=True)
add_stat_annotation(ax, data=box_data, x="type", y="Spike cavity",box_pairs=[("Fab", "Nb")], test='t-test_welch', text_format='star', loc='outside', verbose=2)
ax.set_xlabel('')


plt.savefig('/cs/usr/tomer.cohen13/lab/Epitop_paper/curvature_boxplot.png', dpi=500)
box_data.to_csv("/cs/usr/tomer.cohen13/lab/Epitop_paper/curvature_data.csv", index=False)
exit()



########################################################################################################################
#                                               Nb Ab compare                                                          #
########################################################################################################################


nb_matrix = nb_epitops_df.drop(["Unnamed: 0"], axis=1)
ab_matrix = epitops_nr.drop(["Unnamed: 0", "index"],axis=1)

#################
nb_matrix['ANTIBODY'] = nb_matrix['ANTIBODY'].map(lambda x: x.split('.')[0].upper() if "RBDtr" not in x else x.split('.')[0])
nb_matrix.rename({"ANTIBODY": "PDB", "PDB": "NB"}, axis=1, inplace=True)
ab_matrix.rename({"ANTIBODY": "PDB", "PDB": "AB"}, axis=1, inplace=True)
nb_matrix.to_csv("/cs/usr/tomer.cohen13/lab/Epitop_paper/nb_epitopes_raw_matrix.csv", index=False)
ab_matrix.to_csv("/cs/usr/tomer.cohen13/lab/Epitop_paper/ab_epitopes_raw_matrix.csv", index=False)
exit()
#################

nb_matrix.index = nb_matrix["PDB"]
nb_matrix = nb_matrix.drop(["PDB","ANTIBODY"], axis=1)

ab_matrix.index = ab_matrix["ANTIBODY"]
ab_matrix = ab_matrix.drop(["PDB","ANTIBODY"], axis=1)

our_nbs = ["Nb17", "Nb20","Nb21","Nb34","Nb36","Nb93","Nb95", "Nb105"]
euclidian = []
euclidian_dist = []
jaccard = []
jaccard_score = []
ops = 0
bad_ab = []
for my_nb in our_nbs:
    nb_vector = np.array(nb_matrix.loc[my_nb])
    best_eucl = ""
    best_eucl_dist = np.inf
    best_jaccard = ""
    best_jaccard_score = 0

    for antibody in ab_matrix.index:
        ab_vector = np.array(ab_matrix.loc[antibody])
        e_dist = np.sum((nb_vector - ab_vector)**2)
        if e_dist < best_eucl_dist:
            best_eucl_dist = e_dist
            best_eucl = antibody

        mone = np.sum(nb_vector*ab_vector)
        mechane = np.sum((nb_vector + ab_vector) > 0)
        j_score = mone/mechane
        if j_score > best_jaccard_score:
            best_jaccard_score = j_score
            best_jaccard = antibody
        if np.sum(ab_vector) == 0:
            ops +=1
            bad_ab.append(antibody)

    euclidian.append(best_eucl)
    euclidian_dist.append(best_eucl_dist)

    jaccard.append(best_jaccard)
    jaccard_score.append(best_jaccard_score)

jaccard_ab = [epitops_nr[epitops_nr["ANTIBODY"] == i]["PDB"].iloc[0] for i in jaccard]
euclidian_ab = [epitops_nr[epitops_nr["ANTIBODY"] == i]["PDB"].iloc[0] for i in euclidian]

score_df = pd.DataFrame({"Nb": our_nbs, "PDB-jaccard": jaccard, "Ab-jaccard": jaccard_ab, "jaccard score": jaccard_score, "PDB-euclidean": euclidian, "Ab-euclidean": euclidian_ab, "Euclidean dist": euclidian_dist})
score_df.to_csv("/cs/usr/tomer.cohen13/lab/Epitop_paper/best_match_nb.csv", index=False)


########################################################################################################################
#                                               Nb figure                                                              #
########################################################################################################################

nb_epitops_df["ANTIBODY"] = nb_epitops_df["ANTIBODY"].map(lambda x: x.split(".")[0].upper() if x != "7jvb.pdb" else x.split(".")[0].upper() + " (Nb20)")
nb_epitops_df["ANTIBODY"] = nb_epitops_df["ANTIBODY"].map(lambda x: x if "NB" not in x else "Nb" + x.split("_")[0].split("NB")[1])

nb_epitops_df.index = [i + "@" + j for i ,j in zip(nb_epitops_df["ANTIBODY"], nb_epitops_df["PDB"])]

nb_epitops_df = nb_epitops_df.drop( ["PDB", "ANTIBODY", "Unnamed: 0"], axis=1)
# nb_epitops_df = nb_epitops_df.drop( ["7JWB"], axis=0)

kmeans = KMeans(n_clusters=3)
y = kmeans.fit_predict(nb_epitops_df)
nb_epitops_df["cluster"] = y
nb_epitops_df = nb_epitops_df.sort_values(by="cluster")
nb_epitops_df = nb_epitops_df.drop(['cluster'], axis=1)

# fig1
fig_nb,ax_nb = plt.subplots(1,1,figsize=(26,10))
heat = sns.heatmap(nb_epitops_df,cmap=[plt.get_cmap('coolwarm')(0), plt.get_cmap('coolwarm')(np.inf)],  # Choose a squential colormap
                annot_kws={'fontsize':11},  # Reduce size of label to fit
                linewidth=0.01,  # Add gridlines
                linecolor="#222",# Adjust gridline color
                ax=ax_nb,        # Arrange in subplot
                cbar_kws={"shrink": 0.5, "pad": 0.01})
colorbar = ax_nb.collections[0].colorbar
colorbar.set_ticks([1 / 4, 3 / 4])
colorbar.set_ticklabels(['False', 'True'])
colorbar.set_label("Contact", fontsize=12)
ax_nb.set_title('Nb epitopes', fontsize=20)
ax_nb.set_ylabel('Nb PDB', fontsize=12)
ax_nb.set_xlabel('Position', fontsize=12)

heat.set_yticklabels([i.split("@")[0] for i in nb_epitops_df.index])

# iterate through both the labels and the texts in the heatmap (ax.texts)
for lab in heat.get_yticklabels():
    text = lab.get_text()
    if text in ["Nb17", "Nb34", "Nb36", "Nb21", "Nb20", "Nb93", "Nb95", "Nb105", "7JVB (Nb20)"]: # highlight
        # set the properties of the ticklabel
        lab.set_weight('bold')
        # lab.set_size(20)
        lab.set_color('green')

plt.savefig('/cs/usr/tomer.cohen13/lab/Epitop_paper/label_pdb.png',dpi=500)


fig_nb,ax_nb = plt.subplots(1,1,figsize=(26,10))
heat = sns.heatmap(nb_epitops_df,cmap=[plt.get_cmap('coolwarm')(0), plt.get_cmap('coolwarm')(np.inf)],  # Choose a squential colormap
                annot_kws={'fontsize':11},  # Reduce size of label to fit
                linewidth=0.01,  # Add gridlines
                linecolor="#222",# Adjust gridline color
                ax=ax_nb,        # Arrange in subplot
                cbar_kws={"shrink": 0.5, "pad": 0.01})
colorbar = ax_nb.collections[0].colorbar
colorbar.set_ticks([1 / 4, 3 / 4])
colorbar.set_ticklabels(['False', 'True'])
colorbar.set_label("Contact", fontsize=12)
ax_nb.set_title('Nb epitopes', fontsize=20)
ax_nb.set_ylabel('Nb', fontsize=12)
ax_nb.set_xlabel('Position', fontsize=12)

heat.set_yticklabels([i.split("@")[1] for i in nb_epitops_df.index])

# iterate through both the labels and the texts in the heatmap (ax.texts)
for lab in heat.get_yticklabels():
    text = lab.get_text()
    if text in ["Nb17", "Nb34", "Nb36", "Nb21", "Nb20", "Nb93", "Nb95", "Nb105", "7JVB (Nb20)"]: # highlight
        # set the properties of the ticklabel
        lab.set_weight('bold')
        # lab.set_size(20)
        lab.set_color('green')

plt.savefig('/cs/usr/tomer.cohen13/lab/Epitop_paper/label_ab.png',dpi=500)


########################################################################################################################
#                                               Plot figure                                                            #
########################################################################################################################

fig, ax = plt.subplots(6, 1, figsize=(26,10))

for i,d in enumerate([consurf_scores, mutations, mean_epitop_nr, mean_nb_epitop_nr, ace2_bind, ace2_contact]):
    robust = False
    if i == 1:
        robust = True
    df = pd.DataFrame({"value":d})
    df.index = positions
    df = df.transpose()
    cmap = "coolwarm"
    if i == 5:
        cmap = [plt.get_cmap('coolwarm')(0), plt.get_cmap('coolwarm')(np.inf)]

    sns.heatmap(df,cmap=cmap,  # Choose a squential colormap
                annot_kws={'fontsize':11},  # Reduce size of label to fit
                linewidth=0.01,  # Add gridlines
                linecolor="#222",# Adjust gridline color
                ax=ax[i],        # Arrange in subplot
                robust=robust,
                cbar_kws={"shrink": 1, "pad": 0.01}, yticklabels=False)
    ax[i].patch.set_alpha(0)
    if i == 5:
        colorbar = ax[i].collections[0].colorbar
        colorbar.set_ticks([1 / 4 , 3/4 ])
        colorbar.set_ticklabels(['False', 'True'])
fig.patch.set_facecolor('white')

ax[0].set_title('ConSurf score')
ax[1].set_title('SARS-CoV-2 mutations')
ax[2].set_title('IgG epitopes')
ax[3].set_title('Nb epitopes')
ax[4].set_title('ACE2 binding')
ax[5].set_title('ACE2 contact')


ax[0].set_ylabel('Score')
ax[1].set_ylabel('Count')
ax[2].set_ylabel('Contact')
ax[3].set_ylabel('Contact')
ax[4].set_ylabel('Binding')
ax[5].set_ylabel('Contact')

ax[0].set_xlabel('Position')
ax[1].set_xlabel('Position')
ax[2].set_xlabel('Position')
ax[3].set_xlabel('Position')
ax[4].set_xlabel('Position')
ax[5].set_xlabel('Position')


plt.tight_layout()
plt.savefig('/cs/usr/tomer.cohen13/lab/Epitop_paper/heatmap.png', dpi=500)

save_df = pd.DataFrame({"position": positions, "consurf_score": consurf_scores, "mutations_count": mutations, "ab_epitopes":mean_epitop_nr, "nb_epitopes": mean_nb_epitop_nr ,"ace2_binding": ace2_bind, "ace2_contact": ace2_contact})
save_df.to_csv("/cs/usr/tomer.cohen13/lab/Epitop_paper/heatmap_data.csv", index=False)

exit()
########################################################################################################################
#                                               All spike Figure                                                       #
########################################################################################################################

consurf_df = pd.read_csv("/cs/labs/dina/tomer.cohen13/Epitop_paper/consurf_all.grades", skipfooter=5, skiprows=list(range(13)) + [14], delim_whitespace=True, skipinitialspace=True)

positions = []
consurf_scores = []
begin = False
prev_pos = 0

for score, pos in zip(consurf_df["SCORE"],consurf_df["3LATOM"]):

    if pos == '-' and not begin:
        continue
    else:
        begin = True
        if pos == "-":
            positions.append(prev_pos+1)
        else:
            positions.append(int(pos[3:-2]))
        prev_pos = positions[-1]
        consurf_scores.append(score)
        if positions[-1] == 1147:
            begin = False

mutations_df = pd.read_csv("/cs/usr/tomer.cohen13/lab/Epitop_paper/coronavirus3d-S-mutations.tsv",sep='\t')
only_ms = mutations_df[mutations_df['Annotation Type'].map(lambda x: 'missense_variant' in x)]

mutations = np.zeros(len(positions), dtype=int)
start_pos = positions[0]
end_pos = positions[-1]

for row in only_ms.index:

    already_visited = set()
    row_mut = only_ms["AA change"][row].split(";")
    for mut in row_mut:
        pos = mut[13:].split("-")
        if len(pos) == 1:
            pos = int(''.join(c for c in pos[0] if c.isdigit()))
            if (start_pos <= pos <= end_pos) and (pos not in already_visited):
                mutations[pos-start_pos] += only_ms["Count"][row]
                already_visited.add(pos)

fig, ax = plt.subplots(4, 1, figsize=(30,10))


for i,d in enumerate([consurf_scores, mutations]):
    robust = False
    if i == 1:
        vmax = np.percentile(d, 98)
        vmin = np.percentile(d, 2)
    else:
        vmax= np.max(d)
        vmin = np.min(d)
    df = pd.DataFrame({"value":d[0:len(d)//2]})
    df.index = positions[0:len(d)//2]
    df = df.transpose()
    sns.heatmap(df,cmap="coolwarm",  # Choose a squential colormap
                annot_kws={'fontsize':11},  # Reduce size of label to fit
                linewidth=0.01,  # Add gridlines
                linecolor="#222",# Adjust gridline color
                ax=ax[i],        # Arrange in subplot
                vmin=vmin, vmax=vmax,
                cbar_kws={"shrink": 1, "pad": 0.01}, yticklabels=False)
    ax[i].patch.set_alpha(0)

    df = pd.DataFrame({"value":d[len(d)//2:]})
    df.index = positions[len(d)//2:]
    df = df.transpose()
    sns.heatmap(df,cmap="coolwarm",  # Choose a squential colormap
                annot_kws={'fontsize':11},  # Reduce size of label to fit
                linewidth=0.01,  # Add gridlines
                linecolor="#222",# Adjust gridline color
                ax=ax[i+2],        # Arrange in subplot
                vmin=vmin, vmax=vmax,
                cbar_kws={"shrink": 1, "pad": 0.01}, yticklabels=False)
    ax[i].patch.set_alpha(0)


fig.patch.set_facecolor('white')

ax[0].set_title('ConSurf score')
ax[1].set_title('SARS-CoV-2 mutations')

ax[0].set_ylabel('Score')
ax[1].set_ylabel('Count')

ax[0].set_xlabel('Position')
ax[1].set_xlabel('Position')


ax[2].set_title('ConSurf score')
ax[3].set_title('SARS-CoV-2 mutations')

ax[2].set_ylabel('Score')
ax[3].set_ylabel('Count')

ax[2].set_xlabel('Position')
ax[3].set_xlabel('Position')

plt.tight_layout()
plt.savefig('/cs/usr/tomer.cohen13/lab/Epitop_paper/heatmap_all.png', dpi=500)


save_df = pd.DataFrame({"position": positions, "consurf_score": consurf_scores, "mutations_count": mutations})
save_df.to_csv("/cs/usr/tomer.cohen13/lab/Epitop_paper/heatmap_all_data.csv", index=False)


