import argparse

import pandas as pd

def load_eggnog(input_file, no_isoform=False):

	#load file
	df = pd.read_csv(input_file, sep='\t', comment='#', header=None) 

	#keep genes and GO columns
	df = df.drop(columns=[1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20]) 

	#rename columns
	df = df.rename(columns={0: 'gene', 9: 'GO'})

	if not no_isoform:
		#remove isoform tag (I used protein longest transcript but forgot to remove tags)
		df.gene = df.gene.str[:-2]

	tmp = df['GO'].str.split(',').apply(pd.Series, 1).stack()
	tmp.index = tmp.index.droplevel(-1)
	tmp.name = 'GO'

	#remove the GO column from previous dataframe
	del df['GO']

	#replace it by the GO from temporary dataframe
	df = df.join(tmp)
	df = df.reset_index()

	#Clean up useless added column
	del df['index']

	# Re-order columns to have GO appear first
	df = df[['GO', 'gene']]

	# Remove genes without any GO annotation
	df = df[df.GO != '-']

	return df


def load_onto(input_file):
	d = {}
	with open(input_file, 'r') as infile:
		for i, line in enumerate(infile):
			if i == 0: continue
			line = line.strip().split('\t')
			d[line[0].strip('"')] = [line[1].strip('"')]
	return d


def load_terms(input_file):
	df = pd.read_csv(input_file, sep='\t')
	return df


def format_outputs(go_annot, goterms, onto, out_go_annot, out_terms, domain='BP'):

	biological_process = [i for i in onto if onto[i]==[domain]]
	to_keep = set(go_annot['GO'].unique()).intersection(set(biological_process))

	go_annot = go_annot[go_annot.GO.isin(to_keep)]
	nb_genes_with_go = len(go_annot['gene'].unique())
	print(f'Number of genes annotated with >= 1 GO: {nb_genes_with_go}')

	go_annot.to_csv(out_go_annot, sep='\t', index=False)

	goterms = goterms[goterms.index.isin(to_keep)]
	goterms = goterms.reset_index()
	goterms.to_csv(out_terms, sep='\t',index=False)




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input', type=str, required=True, help='eggnog-mapper file')

    parser.add_argument('--onto', type=str,  required=True)

    parser.add_argument('--terms', type=str,  required=True)

    parser.add_argument('-o', '--out', type=str, required=True)

    parser.add_argument('--out_terms', type=str, required=True)

    parser.add_argument('--no_isoforms', action='store_true')

    parser.add_argument('--domain', type=str, required=False, default='BP')

    args = vars(parser.parse_args())


    egg = load_eggnog(args['input'], args['no_isoforms'])

    onto = load_onto(args['onto'])

    goterms = load_terms(args['terms'])

    format_outputs(egg, goterms, onto, args['out'], args['out_terms'], args['domain'])