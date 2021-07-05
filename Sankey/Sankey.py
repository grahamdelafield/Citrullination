import re
import os
import copy
import pandas as pd
import holoviews as hv
import matplotlib.pyplot as plt
hv.extension("matplotlib")
import warnings
warnings.filterwarnings("ignore")


def read_protein_file(filename, sheet_name, flag=True):
    file = pd.ExcelFile(filename)
    df = pd.DataFrame()
    for sheet in file.sheet_names:
        if sheet=='All IDs':
            sub = pd.read_excel(filename, sheet_name=sheet, usecols=[0]).dropna()
        else:
            sub = pd.read_excel(filename, sheet_name=sheet, usecols=[0]).dropna()
        sub.columns = [sheet]
        if df.empty:
            df = sub
        else:
            df = df.join(sub, how='outer')
    if flag:
        return df


def parse_GO_file(filename, sep='\t'):
    """As many of the modified proteins have not previously been reported as
       citrullinated, we wish to provide analysis of which proteins are
       currently reported as modified and/or citrullinated in Uniprot.

       The file used for the subsequent code was obtained by searching each
       protein accession number in Uniprot and modifying the column results to
       display PTMs."""

    go_df = pd.read_csv(filename, sep)
    cols = [col for col in go_df]
    cols[-1] = "PTM"
    go_df.columns = cols
    go_df = go_df.fillna(0)

    previous_report = []
    go_query = go_df[['Entry', 'PTM']]

    pattern = r'[Cc]itrull'
    for item in go_query.PTM.tolist():
        if item == 0:
            previous_report.append('No Prior PTM in Uniprot')
        # elif re.search(pattern, item):
        #     previous_report.append("Reported as Citrullinated")
        else:
            previous_report.append("Uniprot PTM but not Citrullinated")

    go_query.loc[:, "Previous_report"] = pd.Series(previous_report)
    return go_query


def get_proteins(dataframe):
    """Reads identified proteins from dataframe and returns
       a list.

       attributes: dataframe (type:pd.dataframe)"""

    all_proteins = dataframe['All IDs'].tolist()
    return all_proteins


def parse_regions(dataframe):
    """In order to distinguish identified proteins based on the tissue in
       which they are found, the brain will be treated as a unique class with
       five distinct tissue types.

       attributes: dataframe (type:pd.dataframe)"""

    brain_regions = [col for col in dataframe.iloc[:, 1:6]]
    organs = [col for col in dataframe if col not in brain_regions
              and col != "All IDs"]
    return brain_regions, organs


def create_count_dict(dataframe, brain_regions, organs):
    """For use in the Sankey Plot later, the number of proteins that exist in
       each region must be determined. Though some values will not appear in
       the final figure, the numbers used here are necessary to establish
       correct weights in the plot.

       attributes: dataframe     (type:pd.dataframe)
                   brain_regions (type: list)
                   organs        (type: list)"""
    all_proteins = get_proteins(dataframe)
    count_d = {}

    for region in brain_regions:
        count = 0
        region_proteins = dataframe[region].tolist()
        for protein in region_proteins:
            if protein in all_proteins:
                count += 1
        count_d[region] = count
    for organ in organs:
        count = 0
        organ_proteins = dataframe[organ].tolist()
        for protein in organ_proteins:
            if protein in all_proteins:
                count += 1
        count_d[organ] = count

    return count_d


def orig_sankey(protein_df, brain_regions, organs, count_dict, go_df):
    """In order to use the Holoviews.Sankey plot, all data must be called into
       dataframe. As stated before, this dataframe will produce numbers that are
       not used in the final figure. Though these values are correct and true,
       they do not account for overlapping occurrences. The values used in this
       dataframe are employed to establish correct weights in the Sankey Plot.

       The true values can be displayed by calling a "print" of the
       corrected_values function.


       attributes: protein_df (type: dataframe)
                   brain_regions (type: list)
                   organs (type: list)
                   count_dict (type: dictionary)
                   go_df (type: dataframe)"""
    s = []
    t = []
    value = []
    brain_total = 0
    body_total = 0
    for key in count_dict:
        if key not in brain_regions:
            body_total += count_dict[key]
            s.append("Body")
            t.append(key)
            value.append(count_dict[key])
        else:
            brain_total += count_dict[key]
            s.append("Brain")
            t.append(key)
            value.append(count_dict[key])
    s.append("Total")
    t.append("Brain")
    value.append(brain_total)
    s.append("Total")
    t.append("Body")
    value.append(body_total)

    for region in brain_regions:
        region_proteins = protein_df[region].tolist()
        sub_df = go_df[go_df.Entry.isin(region_proteins)]
        val_counts = sub_df.Previous_report.value_counts()
        for item in val_counts.items():
            s.append(region)
            t.append(item[0])
            value.append(item[1])

    for organ in organs:
        organ_proteins = protein_df[organ].tolist()
        sub_df = go_df[go_df.Entry.isin(organ_proteins)]
        val_counts = sub_df.Previous_report.value_counts()
        for item in val_counts.items():
            s.append(organ)
            t.append(item[0])
            value.append(item[1])

    d = {}
    d['source'] = pd.Series(s)
    d['target'] = pd.Series(t)
    d['value'] = pd.Series(value)
    return pd.DataFrame(d)

def prop_sankey(brain_regions, organs, count, corrected, protein_df, go_df):
    s, t = [], []
    values = []
    # put in total, brain and organs
    for group in ['Brain', 'Body']:
        s.append('Total')
        t.append(group)
        values.append(corrected[group])
    for r in brain_regions:
        s.append('Brain')
        t.append(r)
        values.append(count[r])
    for r in organs:
        s.append('Body')
        t.append(r)
        values.append(count[r])
    # put in GO
    for region in brain_regions:
        region_proteins = protein_df[region].tolist()
        sub_df = go_df[go_df.Entry.isin(region_proteins)]
        val_counts = sub_df.Previous_report.value_counts()
        for item in val_counts.items():
            s.append(region)
            t.append(item[0])
            values.append(item[1])

    for organ in organs:
        organ_proteins = protein_df[organ].tolist()
        sub_df = go_df[go_df.Entry.isin(organ_proteins)]
        val_counts = sub_df.Previous_report.value_counts()
        for item in val_counts.items():
            s.append(organ)
            t.append(item[0])
            values.append(item[1])

    d = {}
    d['source'] = pd.Series(s)
    d['target'] = pd.Series(t)
    d['value'] = pd.Series(values)
    return pd.DataFrame(d)

def plot_sankey(sankey_dataframe, filename, save=True):
    """The plotting code below is credited to the Holoviews user Gallery:
    http://holoviews.org/gallery/demos/bokeh/energy_sankey.html#bokeh-gallery-energy-sankey

    By default this function will save a PNG of the plot used. If not desired,
    set the "save" kwarg to False or change the hv.save() extension to the
    desired format."""

    sankey = hv.Sankey(sankey_dataframe, label='Citrulinated Proteins')
    sankey.opts(label_position='right', edge_color='target', node_color='index', cmap='tab20c')
    if save:
        # hv.save(sankey, filename+'.png')
        hv.save(sankey, filename+'.svg')

def unique(df, columns):
    g =  ['Brain', 'Body', 'Total']
    l = [5, 6, 11]
    res = {}
    all_c = []
    for c in columns:
        all_c.extend(c)
    columns.append(all_c)
    output = []
    for group in columns:
        print(group)
        k = g[l.index(len(group))]
        sub = df.loc[:, group]
        sub = sub.melt()
        sub = sub.dropna()
        total = sub.value.unique()
        res[k] = len(total)
        text = f'Group: "{group}" has {len(total)} total proteins'
        output.append(text)
    with open('Sankey Corrected Values.txt', 'a') as f:
        f.write('\n'.join(output))
        f.write('\n')

    return res



def corrected_values(go_df, protein_df, unique_dict):
    sliced_go = go_df[go_df.Entry.isin(protein_df['All IDs'])]
    print('sliced', len(sliced_go))
    d = {}
    for item in sliced_go.Previous_report.tolist():
        if item not in d:
            d[item] = 1
        else:
            d[item] += 1
    output = []
    for k, v in d.items():
        unique_dict[k] = v
        text = f'Group: "{k}" has {v} total proteins'
        output.append(text)
    with open('Sankey Corrected Values.txt', 'a') as f:
        f.write('\n'.join(output))
    return unique_dict


def plot_site_overlap(filename, sheet_name):
    df = pd.read_excel('Sites with different modifications.xlsx',
                       'Citrullination sites')

    targets = [col for col in df] + ['No Overlap']

    all_overlap = []
    counts = []
    cit_d = {}

    for item in df.iloc[:, 0]:
        cit_d[item] = None

    for i in df.iloc[:, 1:]:
        count = 0
        subset = set(df[i].tolist())
        for site in subset:
            if site in cit_d:
                count += 1
                all_overlap.append(site)
        counts.append(count)
    count = 0
    for item in all_overlap:
        if item not in cit_d:
            count += 1
    counts.append(len(cit_d)-len(all_overlap))
    names = targets[1:]

    my_circle = plt.Circle((0, 0), 0.7, color='white')

    # print(counts)

    plt.rcParams['figure.figsize'] = (10, 10)
    plt.pie(counts, labels=names, colors=['red', 'green', 'blue', 'skyblue'])
    p = plt.gcf()
    p.gca().add_artist(my_circle)
    p.savefig('Donut.png')


def main():
    protein_df = read_protein_file('../Tissue-specific Cit IDs.xlsx', 'All IDs')
    print(protein_df)
    brain_regions, organs = parse_regions(protein_df)
    print(brain_regions, organs)
    ud = unique(protein_df, [brain_regions, organs])
    print('UD', ud)
    go_query = parse_GO_file('./Citrullinated_Uniprot.tab')
    # go_query.to_csv('GO.csv')
    count_dict = create_count_dict(protein_df, brain_regions, organs)
    print(count_dict)
    corrected = corrected_values(go_query, protein_df, ud)
    print(corrected)
    prop_sank = prop_sankey(brain_regions, organs, count_dict, corrected, protein_df, go_query)
    or_sank = orig_sankey(protein_df, brain_regions, organs,
                                   count_dict, go_query)

    plot_sankey(prop_sank, 'Proportional')
    plot_sankey(or_sank, 'Original')

if __name__ == '__main__':
    print(os.getcwd())
    main()