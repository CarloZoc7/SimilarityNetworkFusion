import gzip
import pandas as pd
import re
import os
import platform

from collections import defaultdict
from tqdm import tqdm
# source: https://gist.github.com/slowkow/8101481

GTF_HEADER  = ['seqname', 'source', 'feature', 'start', 'end', 'score',
			   'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')

"""
	Class to map the ENG.
	The source function to read and manage the .GTF file is taken from the following repository:
		https://gist.github.com/slowkow/8101481
"""

class MappingENG:

	def dataframe_gtf(self, filename):
	# Open an optionally gzipped GTF file and return a pandas.DataFrame.
	# Each column is a list stored as a value in this dict.

		result = defaultdict(list)

		for i, line in enumerate(self.lines(filename)):
			for key in line.keys():
				# This key has not been seen yet, so set it to None for all
				# previous lines.
				if key not in result:
					result[key] = [None] * i

			# Ensure this row has some value for each column.
			for key in result.keys():
				result[key].append(line.get(key, None))

		return pd.DataFrame(result)

	def mapENG_GEN(self, genes):
		dict_eng_gen = dict()
		for gene in tqdm(genes):
			if '.' in gene:
				gene_id = gene.split('.')[0]
				gene_version = gene.split('.')[1]

				result = self.df_gtf[self.df_gtf['gene_id'] == gene_id]
				result = result[result['gene_version'] == gene_version]
				if result.empty is True:
					 # looking for gene_id without version
					result = self.df_gtf[self.df_gtf['gene_id'] == gene_id]
			else:
				 result = self.df_gtf[self.df_gtf['gene_id'] == gene]
			
			gene_name = gene # of default I store the eng name
			if result.empty is False:
				gene_name = result['gene_name'].iloc[0] #take the first one
			dict_eng_gen[gene] = gene_name
		return dict_eng_gen

	def fit(self, X, y=None):
		return self

	def transform(self, X, y=None):
		# understand if the no mapping features have to drop

		if os.path.exists(self.path_to_mappedPKL) is False:
			mydict = self.mapENG_GEN(X.columns)
			del self.df_gtf
			X.rename(columns=mydict, errors="raise", inplace=True)
			X.to_pickle(self.path_to_mappedPKL)
			return X
		else:
			return self.X

	def fit_transform(self, X, y=None):
		return self.fit(X, y).transform(X, y)

	def lines(self, filename):
		# Open an optionally gzipped GTF file and generate a dict for each line.
		fn_open = gzip.open if filename.endswith('.gz') else open

		with fn_open(filename) as fh:
			for line in fh:
				if line.startswith('#'):
					continue
				else:
					yield self.parse(line)


	def parse(self, line):
		# Parse a single GTF line and return a dict.
		result = {}

		fields = line.rstrip().split('\t')

		for i, col in enumerate(GTF_HEADER):
			result[col] = self._get_value(fields[i])

		# INFO field consists of "key1=value;key2=value;...".
		infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

		for i, info in enumerate(infos, 1):
			# It should be key="value".
			try:
				key, _, value = re.split(R_KEYVALUE, info, 1)
			# But sometimes it is just "value".
			except ValueError:
				key = 'INFO{}'.format(i)
				value = info
			# Ignore the field if there is no value.
			if value:
				result[key] = self._get_value(value)

		return result


	def _get_value(self, value):
		if not value:
			return None

		# Strip double and single quotes.
		value = value.strip('"\'')

		# Return a list if the value has a comma.
		if ',' in value:
			value = re.split(R_COMMA, value)
		# These values are equivalent to None.
		elif value in ['', '.', 'NA']:
			return None

		return value

	def __init__(self, name=''):
		self.path_to_mappedPKL = './data-ready/data_mapped_'+ name +'.pkl'
		self.path_to_dfPKL = './data-ready/df_gtf.pkl'
		self.path_to_homoGTF = './data-ready/Homo_sapiens.GRCh37.87.gtf'

		if platform.system() == 'Windows':
			self.path_to_mappedPKL = self.path_to_mappedPKL.replace('/', '\\')
			self.path_to_dfPKL = self.path_to_dfPKL.replace('/', '\\')
			self.path_to_homoGTF = self.path_to_homoGTF.replace('/', '\\')

		if os.path.exists(self.path_to_mappedPKL) is True:
			self.X = pd.read_pickle(self.path_to_mappedPKL)
		else:
			if os.path.exists(self.path_to_dfPKL) is True: 
				self.df_gtf = pd.read_pickle(self.path_to_dfPKL)
			else:
				self.df_gtf = self.dataframe_gtf(self.path_to_homoGTF)
				self.df_gtf.to_pickle(self.path_to_dfPKL)


