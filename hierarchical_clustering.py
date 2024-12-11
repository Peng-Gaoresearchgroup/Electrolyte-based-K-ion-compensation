import os
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, AllChem, Descriptors3D, Lipinski
from rdkit.ML.Cluster import Butina
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Draw
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Descriptors3D import Asphericity, Eccentricity, InertialShapeFactor, NPR1, NPR2, PMI1, PMI2, PMI3,RadiusOfGyration, SpherocityIndex
from sklearn import preprocessing
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram,fcluster,linkage
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.pyplot import cm
import matplotlib as mpl
from collections import defaultdict

SEED=42
# AL == Anode Limit, SE == Solubility Energy.
# 'get','2d','3d' store the descriptors names. 
# 'distance_threshold' is a hyperparameters for AgglomerativeClustering.
# 'color_threshold' defines the loacation of partition line in output tif images. It helps in the subsequent presentation of the data without computational nature. 
setting_AL={'get':['NumAtoms','NumHeavyAtoms'],
                '2d':['NumHeteroatoms','NumRotatableBonds','RingCount','NumAromaticRings','NumAromaticCarbocycles','MolLogP','MolMR'],
                '3d':['Asphericity','Eccentricity','NPR1','NPR2','PMI1','PMI2','PMI3','RadiusOfGyration','SpherocityIndex'],
                'distance_threshold':0.5,
                'color_threshold':1.25}
setting_SE={'get':['NumAtoms','NumHeavyAtoms'],
                '2d':['NumHeteroatoms','NumRotatableBonds','RingCount','NumAromaticRings','NumAromaticCarbocycles','NumHAcceptors','MolWt','HeavyAtomMolWt','TPSA'],
                '3d':['Asphericity','Eccentricity','NPR1','NPR2','PMI1','PMI2','PMI3','RadiusOfGyration','SpherocityIndex'],
                'distance_threshold':0.5,
                'color_threshold':1.88}
jobs={'AL':setting_AL,'SE':setting_SE}
fig_size,dpi,font,font_size,boder_linewidth,tree_linewidth=(10,2),400,'Arial',5,0.25,0.5

def get_desdescriptors(smiles,job,seed=SEED):
    '''Get descriptors from a SMILES string.'''
    mol=Chem.MolFromSmiles(smiles) 
    if mol is None:
        raise ValueError("SMILES Error")
    # Choose a job
    descriptors_dic=jobs[job]
    # Get 2D descriptors
    descriptors_value={f'{i}': getattr(mol, f'Get{i}')() for i in descriptors_dic['get']}
    descriptors_value.update({f'{i}': getattr(Descriptors, f'{i}')(mol) for i in descriptors_dic['2d']})
    # 3D Embed
    AllChem.EmbedMolecule(mol, useRandomCoords=True,maxAttempts=5000, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(mol)
    mol=Chem.AddHs(mol)
    params = AllChem.ETKDG()
    params.randomSeed = seed
    AllChem.EmbedMolecule(mol, params)
    # Get 3D descriptors
    descriptors_value.update({f'{i}': getattr(Descriptors3D, f'{i}')(mol) for i in descriptors_dic['3d']})
    return pd.Series(descriptors_value)

def get_dendrogram(model, p, line_thickness, num_clusters=4,colour_threshold=0):
    '''Get dendrogram information(but not draw it) form an AgglomerativeClustering instance'''
    class_dict = {}
    counts= []
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1
            else:
                current_count += counts[child_idx - n_samples]
        counts.append(current_count)
    
    Z = np.column_stack([model.children_, model.distances_, counts]).astype(float)
    hierarchy.set_link_color_palette(['#DC143C', '#FF8C00', '#2E8B57', '#1E90FF', 'r'])
    dendrogram_results = dendrogram(Z, no_plot=True, p=p, truncate_mode='lastp', show_contracted=True,color_threshold=colour_threshold, above_threshold_color='grey')
    classes = fcluster(Z, t=num_clusters, criterion='distance')
    class_dict = defaultdict(list)
    for index, cls in enumerate(classes):
        class_dict[cls].append(index)
    for cls, idx_list in class_dict.items():
        print(f"Class {cls}: Indices {idx_list}")
    for x, y, color in zip(dendrogram_results['icoord'], dendrogram_results['dcoord'], dendrogram_results['color_list']):
        plt.plot(x, y, lw=line_thickness, color=color)
    plt.ylim(0, Z[:, 2].max() * 1.05)
    plt.xlim(0, 10 * len(dendrogram_results['ivl']))
    return dendrogram_results["ivl"]

def plt_cluster(model,data_len,fig_save,fig_size,dpi,linewidth,font_size,font,color_threshold,tree_linewidth):
    '''Draw dendrogram form an AgglomerativeClustering instance'''
    plt.figure(figsize=fig_size)
    trees = get_dendrogram(model=model, p=data_len,line_thickness=tree_linewidth,colour_threshold=color_threshold)
    plt.ylabel('Distance', fontsize=font_size, fontname=font, weight='bold')
    plt.tick_params(axis='x', labelsize=font_size, width=linewidth, bottom=True, labelbottom=True)
    plt.tick_params(axis='y', width=linewidth, labelsize=font_size)
    plt.xticks(ticks=[5 + 10 * i for i in range(len(trees))],
               labels=[int(node) for node in trees],
               fontsize=font_size, rotation=90)

    ax = plt.gca()
    ax.yaxis.set_ticks_position('right')
    ax.yaxis.set_label_position('right')
    ax.tick_params(axis='y', labelsize=font_size, labelrotation=90)
    
    for spine in ['top', 'right', 'bottom', 'left']:
        ax.spines[spine].set_linewidth(linewidth)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    plt.axhline(y=color_threshold, color='black', linestyle='--', linewidth=tree_linewidth)
    
    plt.tight_layout()
    os.makedirs(os.path.dirname(fig_save), exist_ok=True)
    plt.savefig(fig_save, format='tif', dpi=dpi, bbox_inches='tight')
    return None

if __name__ == '__main__': 
    for job in jobs:
        data = pd.read_excel('./Data/molecules.xlsx')
        data['SMILES'] = data['SMILES'].apply(lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x)))
        data = pd.concat([data, data['SMILES'].apply(lambda x: get_desdescriptors(smiles=x, job=f'{job}'))], axis=1)
        
        #Data Normalization
        scaler = MinMaxScaler()
        scaled_data = pd.DataFrame(scaler.fit_transform(data.iloc[:, 1:]), columns=data.columns[1:])

        data = pd.concat([data.iloc[:, :1], scaled_data], axis=1)
        data.to_excel(f'./Output/data_{job}.xlsx', index=False)
        model = AgglomerativeClustering(distance_threshold=jobs[job]['distance_threshold'], n_clusters=None).fit(data.iloc[:, 1:])
        plt_cluster(model, data_len=len(data.iloc[:, 1:]), fig_save=f'./Output/cluster_{job}.tif',fig_size=fig_size,dpi=dpi,linewidth=boder_linewidth,font_size=font_size,font=font,color_threshold=jobs[job]['color_threshold'],tree_linewidth=tree_linewidth)
