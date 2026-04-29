####################################################################
###These functions were made possible thanks to Dr. Victor Perez###
####################################################################

import os
import pandas as pd
from pathlib import Path
import numpy as np


#fix the macsiq format

def concat_MACSiq(work_dir,preprocessed_dir, lineage_info):
    files = [x for x in work_dir.iterdir() if x.is_file() and x.suffix == '.csv' and x.name != '1 all_info.csv']
    output_dir = Path(preprocessed_dir/work_dir.name)
    if os.path.exists(output_dir):
        pass 
    else:
        output_dir.mkdir(parents=True, exist_ok=True)
    output_file = Path(output_dir / 'macsiq_phenotyping.csv')
    if os.path.exists(output_file):
        print("File exists")
        return pd.read_csv(output_file)
    else:
        pass
    data=[]
    for f in files:
        df=pd.read_csv(f)
        symbol='@'+df['Cell Id'][0].split('@')[1]
        df.insert(1,'cell_type',df.shape[0]*[f.name[0:-4]])
        df.insert(1,'cell_id',''.join(df['Cell Id'].values).split(symbol)[0:-1])
        data.append(df)
    
    df=pd.concat(data,ignore_index=True)
    df["cell_id"] = pd.to_numeric(df["cell_id"])
    df.sort_values(by=['cell_id'],inplace=True,ignore_index=True)
    df_macsiq=df.loc[:,['cell_id','cell_type']]
    df_macsiq['cell_type'] = df_macsiq['cell_type'].str.replace(r'^\d+\s+', '', regex=True)

    level=[]
    source=[]
    cell_type_label=[]
    for c in df_macsiq.cell_type.values:
        info=lineage_info[c]
        level.append(info[0])
        source.append(info[1])
        cell_type_label.append(info[2])

    df_lineage=pd.DataFrame({'level':level,'parent_cell':source,'type_no': cell_type_label })
    df_macsiq=pd.concat([df_macsiq,df_lineage],ignore_index=False,axis=1)

    df_macsiq.to_csv(output_file,index=False)

    return df_macsiq

##for final phenotyping labels adjustment

def filter_MACSiq(work_dir, preprocessed_dir, pheno_depth):
    output_dir = Path(preprocessed_dir/work_dir.name)
    input_file = Path(output_dir / 'macsiq_phenotyping.csv')
    output_file = Path(output_dir/'final_phenotyping_labels.csv')
    if os.path.exists(output_file):
        print(f"{output_file.stem} already exists")
        return pd.read_csv(output_file)
    else:
        pass
    print(f"Using {input_file.stem} to create {output_file.stem}")
    df_macsiq = pd.read_csv(input_file)
    df_macsiq_filt = df_macsiq.loc[df_macsiq['level']<=pheno_depth]
    repeated_labels=[cellID for cellID, rep in (df_macsiq_filt['cell_id'].value_counts()>1).items() if rep==True]
    repeated_labels.sort()
    remove_indices=[]
    for cellID in repeated_labels:
        cell_lineage=df_macsiq_filt.loc[df_macsiq_filt.cell_id==cellID].level
        index=np.argmax(cell_lineage.values)
        final_type_index=cell_lineage.index.values[index]
        remove_elements=np.setdiff1d(cell_lineage.index.values,[final_type_index]).tolist()
        if remove_elements:
            remove_indices.extend(remove_elements)
    df_macsiq_filt=df_macsiq_filt.drop(index=remove_indices)
    df_macsiq_filt.to_csv(output_file,index=False)
    return df_macsiq_filt

##for Metadata

def meta_MACSiq(work_dir, preprocessed_dir):
    output_dir = Path(preprocessed_dir/work_dir.name)
    output_file = Path(output_dir / 'metadata.csv')
    if os.path.exists(output_file):
        print(f"{output_file.stem} already exists")
        return pd.read_csv(output_file)
    else:
        pass
    print(f"Creating metadata file {output_file.stem}") 
    #load the data with the metadata columns only (79 columns)
    meta = pd.read_csv(work_dir / '1 all_info.csv', usecols=range(1,79))
    #add cell_id column
    meta.insert(0, "cell_id", range(1, 1 + len(meta)))
    #loading the phenotyping data
    df_macsiq_filt = pd.read_csv(output_dir/'final_phenotyping_labels.csv')
    #merge the metadata with the phenotyping data
    meta_joined = df_macsiq_filt.merge(meta, how="left")
    #export the metadata file
    meta_joined.to_csv(output_file,index=False)
    return meta_joined

###for expression
def exp_MACSiq(work_dir, metadata_dir, preprocessed_dir):
    output_dir = Path(preprocessed_dir/work_dir.name)
    output_file = Path(output_dir / 'exp.csv')
    run = work_dir.name.split('_')[0]
    if os.path.exists(output_file):
        print(f"{output_file.stem} already exists")
        return pd.read_csv(output_file)
    else:
        pass
    print("Creating expression file", output_file)
    #load the necessary antibody list
    ab_exp = pd.read_csv(metadata_dir/f"{run}_Ab_panel.csv")
    #load the expression data
    all_columns = pd.read_csv(work_dir / '1 all_info.csv', nrows=0).columns
    selected_columns = all_columns[all_columns.str.contains("Cell Exp")]
    #Columns must not be DAPI, APC, PE 
    selected_columns = [col for col in selected_columns if not any(x in col for x in ['DAPI', 'APC', 'PE', 'FITC'])]
    exp = pd.read_csv(work_dir / '1 all_info.csv', usecols=selected_columns)
    #remove ' Cell Exp' from column names
    exp.columns = exp.columns.str.replace(' Cell Exp', "")
    #Filter out the columns for only ab_exp columns
    available_cols = [col for col in exp.columns if col in ab_exp["exp_col"].values]
    ab_exp = ab_exp[ab_exp["exp_col"].isin(available_cols)]
    exp = exp.loc[:,available_cols]
    #reorder ab_exp
    ab_exp = ab_exp.sort_values(by="exp_col", key=lambda x: x.map({col: i for i, col in enumerate(available_cols)}))
    #Fix the names of the antibodies
    exp.columns = ab_exp["antibody_name"].values
    #Filter out the columns with low staining quality
    low_quality = ab_exp[ab_exp["staining quality"] == "LOW"]["antibody_name"]
    exp = exp.drop(columns=low_quality, errors='ignore')
    #drop low quality from ab_exp
    ab_exp = ab_exp[~ab_exp["antibody_name"].isin(low_quality)]
    #add cell_id column
    exp.insert(0, 'cell_id', range(1, 1 + len(exp)))
    #Only keep rows with a cell type (Remove trash cells)
    meta = pd.read_csv(output_dir / 'metadata.csv')
    cell_ids = meta['cell_id']
    exp = exp[exp['cell_id'].isin(cell_ids)]
    #resort the columns
    exp = exp[['cell_id'] + natsorted(exp.columns[1:])]
    #export the expression data
    exp.to_csv(output_file,index=False)
    return exp
