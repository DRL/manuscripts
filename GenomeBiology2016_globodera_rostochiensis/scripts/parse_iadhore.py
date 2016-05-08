#! /usr/bin/env python
from __future__ import division
import sys
import random
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib as mat
from matplotlib import cm
from matplotlib.ticker import NullFormatter
from matplotlib.lines import Line2D
from matplotlib.colors import rgb2hex
mat.use('agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches

'''
0. parse bed, get protein-contig
1. parse anchorpoints.txt, get contigs, get whether was twisted, account for genes in scaffolds
2. print contig gros contig gpal count genes, twisted?
'''

class Multiplicon():
	def __init__(self, multiplicon):
		self.name = multiplicon
		self.counts = {'baseclusters' : 0, 'anchorpoints' : 0 }
		self.baseclusters = {}
		self.gros_scaffold = ''
		self.gpal_scaffold = ''

def parseBed(bedfile):
	scaffold_of_gene = {}
	genes_of_scaffold = {}
	regions_of_scaffolds = {}
	prev_end = 0
	prev_scaffold = ''
	gene_idx = 1
	with open(bedfile) as fh:
		for l in fh:
			temp = l.rstrip("\n").split("\t")
			if temp[7] == 'gene':
				scaffold = temp[0]
				gene = temp[3]
				gene_start = int(temp[1])
				gene_end = int(temp[2])
				scaffold_of_gene[gene] = (scaffold, gene_idx)
				genes_of_scaffold[scaffold] = genes_of_scaffold.get(scaffold, 0) + 1
				if not scaffold in regions_of_scaffolds: # see scaffold for the first time
					gene_idx = 1
					regions_of_scaffolds[scaffold] = {'genes' : []}
					if (prev_scaffold): # if there is a previous scaffold, deal with it
						prev_region_seq = assembly_d[prev_scaffold][prev_end:len(assembly_d[prev_scaffold])]
						#regions_of_scaffolds[prev_scaffold]['genes'].append(('intergenic', prev_region_seq.count('N'), prev_end, len(assembly_d[prev_scaffold])))
						regions_of_scaffolds[prev_scaffold]['genes'].append(str(prev_region_seq.count('N')))
					if 0 < gene_start: # if gene does not start at beginning of scaffold
						prev_region_seq = assembly_d[scaffold][0:gene_start]
						#regions_of_scaffolds[scaffold]['genes'].append(('intergenic', prev_region_seq.count('N'), 0, gene_start))
						regions_of_scaffolds[scaffold]['genes'].append(str(prev_region_seq.count('N')))
				else: # see scaffold for the second time or more
					prev_region_seq = assembly_d[scaffold][prev_end:gene_start]
					#regions_of_scaffolds[scaffold]['genes'].append(('intergenic', prev_region_seq.count('N'), prev_end, gene_start))
					regions_of_scaffolds[scaffold]['genes'].append(str(prev_region_seq.count('N')))
				gene_region_seq = assembly_d[scaffold][gene_start:gene_end]
				#regions_of_scaffolds[scaffold]['genes'].append((gene, gene_region_seq.count('N'), gene_start, gene_end))
				regions_of_scaffolds[scaffold]['genes'].append(gene)
				prev_end = gene_end
				prev_scaffold = scaffold
				gene_idx += 1

	return scaffold_of_gene, genes_of_scaffold, regions_of_scaffolds

def parseAnchorpoints(anchorpoints_f):
	data = {}
	gros_gene = ''
	gpal_gene = ''
	with open(anchorpoints_f) as fh:
		for l in fh:
			temp = l.rstrip("\n").split("\t")
			if not temp[0] == "id":
				multiplicon = temp[1]
				basecluster = temp[2]
				gene_x = temp[3]
				gene_y = temp[4]
				if not multiplicon in data:
					data[multiplicon] = Multiplicon(multiplicon)
				if not basecluster in data[multiplicon].baseclusters:
					data[multiplicon].baseclusters[basecluster] = {}
					data[multiplicon].counts['baseclusters'] += 1
				if gene_x in scaffold_of_gene_gros:
					gros_gene = gene_x
					gpal_gene = gene_y
				else:
					gros_gene = gene_y
					gpal_gene = gene_x
				gros_scaffold, gros_gene_idx = scaffold_of_gene_gros[gros_gene]
				gpal_scaffold, gpal_gene_idx = scaffold_of_gene_gpal[gpal_gene]
				data[multiplicon].gros_scaffold = gros_scaffold
				data[multiplicon].gpal_scaffold = gpal_scaffold
				data[multiplicon].counts['anchorpoints'] += 1

				if not gros_scaffold in regions_of_scaffolds_d[gpal_scaffold]:
					regions_of_scaffolds_d[gpal_scaffold][gros_scaffold] = ['O' for elem in regions_of_scaffolds_d[gpal_scaffold]['genes']]
				for idx, block in enumerate(regions_of_scaffolds_d[gpal_scaffold]['genes']):
					if block == gpal_gene:
						temp_scaffold, gene_idx = scaffold_of_gene_d[gros_gene]
						regions_of_scaffolds_d[gpal_scaffold][gros_scaffold][idx] = '%s/%s' % (gene_idx, genes_of_scaffold_d[temp_scaffold])

				if not gpal_scaffold in regions_of_scaffolds_d[gros_scaffold]:
					regions_of_scaffolds_d[gros_scaffold][gpal_scaffold] = ['O' for elem in regions_of_scaffolds_d[gros_scaffold]['genes']]
				for idx, block in enumerate(regions_of_scaffolds_d[gros_scaffold]['genes']):
					if block == gros_gene:
						temp_scaffold, gene_idx = scaffold_of_gene_d[gpal_gene]
						regions_of_scaffolds_d[gros_scaffold][gpal_scaffold][idx] = '%s/%s' % (gene_idx, genes_of_scaffold_d[temp_scaffold])

	#for scaffold in regions_of_scaffolds_d:
	#	print scaffold
	#	print 'genes', regions_of_scaffolds_d[scaffold]['genes']
	#	for key in regions_of_scaffolds_d[scaffold]:
	#		if not key == 'genes':
	#			print key, regions_of_scaffolds_d[scaffold][key]
	return data

def parseBaseclusters(baseclusters_f):

	with open("parse_iadhore.blocks.txt", 'w') as fh_w:
		fh_w.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('multiplicon', 'number_of_anchorpoints', 'gros_scaffold', 'gpal_scaffold', 'baseclusters', 'anchorpoints', 'twisted'))
		with open(baseclusters_f) as fh:
			for l in fh:
				temp = l.rstrip("\n").split("\t")
				if not temp[0] == "id":
					multiplicon = temp[1]
					number_of_anchorpoints = temp[2]
					was_twisted = temp[4]
					if was_twisted == '0':
						fh_w.write("%s\t%s\t%s\t%s\t%s\t%s\tno\n" % (multiplicon, number_of_anchorpoints, data[multiplicon].gros_scaffold, data[multiplicon].gpal_scaffold, data[multiplicon].counts['baseclusters'], data[multiplicon].counts['anchorpoints'] ))
					else:
						fh_w.write("%s\t%s\t%s\t%s\t%s\t%s\tyes\n" % (multiplicon, number_of_anchorpoints, data[multiplicon].gros_scaffold, data[multiplicon].gpal_scaffold, data[multiplicon].counts['baseclusters'], data[multiplicon].counts['anchorpoints']))

def getConnections(scaffold, start_scaffold):
	for connection in dict_of_connections[scaffold]:
		cons_of_connection = dict_of_connections[connection]
		if len(cons_of_connection) == 1 and scaffold in cons_of_connection:
			print "%s\t<=>\t%s" % (scaffold, connection)
		else:
			for con_of_cons in cons_of_connection:
				if start_scaffold == con_of_cons:
					print "%s\t<->\t%s" % (start_scaffold, con_of_cons)
				else:
					getConnections(con_of_cons, start_scaffold)
def parseLength(length_f):
	length_dict = {}
	with open(length_f) as fh:
		for l in fh:
			temp = l.rstrip("\n").split("\t")
			scaffold = temp[0]
			length = int(temp[1])
			length_dict[scaffold] = length
	return length_dict

def parseFasta(infile):
	fasta_d = {}
	with open(infile) as fh:
		header, seqs = '', []
		for l in fh:
			if l[0] == '>':
				if (header):
					fasta_d[header] = ''.join(seqs)
				header, seqs = l.lstrip(">").rstrip("\n"), [] # Header is split at first whitespace
			else:
				seqs.append(l.rstrip("\n"))
		fasta_d[header] = ''.join(seqs)
	return fasta_d

class Cluster():
	def __init__(self, idx):
		self.idx = idx
		self.gros_nodes = []
		self.gpal_nodes = []
		self.gros_lengths = []
		self.gpal_lengths = []

if __name__ == "__main__":
	try:
		gros_bed_f = sys.argv[1]
		gros_length_f = sys.argv[2]
		gpal_bed_f = sys.argv[3]
		gpal_length_f = sys.argv[4]
		anchorpoints_f = sys.argv[5]
		baseclusters_f = sys.argv[6]
		gros_assembly_f = sys.argv[7]
		gpal_assembly_f = sys.argv[8]
	except:
		sys.exit("Usage: ./parse_iadhore.py [GROS_BED] [GROS_LENGTH] [GPAL_BED] [GPAL_LENGTH] [ANCHORPOINTS] [BASECLUSTERS] [GROS_ASSEMBLY] [GPAL_ASSEMBLY] \n\n\
	[GROS_BED] 	: Bed file of G. rostochiensis\n\
	[GROS_LENGTH] 	: Length file of G. rostochiensis\n\
	[GPAL_BED] 	: Bed file of G. pallida\n\
	[GPAL_LENGTH] 	: Length file of G. pallida\n\
	[ANCHORPOINTS] 	: Anchorpoints file provided by i-ADHoRe\n\
	[BASECLUSTERS] 	: Baseclusters file provided by i-ADHoRe\n\
	[GROS_ASSEMBLY] : Assembly file of G. rostochiensis\n\
	[GPAL_ASSEMBLY] : Assembly file of G. pallida\n")

	length_of_gros = parseLength(gros_length_f)
	length_of_gpal = parseLength(gpal_length_f)
	lengths = dict(length_of_gpal.items() + length_of_gros.items())
	gros_assembly_d = parseFasta(gros_assembly_f)
	gpal_assembly_d = parseFasta(gpal_assembly_f)
	assembly_d = dict(gros_assembly_d.items() + gpal_assembly_d.items())
	scaffold_of_gene_gros, genes_of_scaffold_gros, regions_of_scaffolds_gros = parseBed(gros_bed_f)
	scaffold_of_gene_gpal, genes_of_scaffold_gpal, regions_of_scaffolds_gpal = parseBed(gpal_bed_f)
	scaffold_of_gene_d = dict(scaffold_of_gene_gpal.items() + scaffold_of_gene_gros.items())
	genes_of_scaffold_d = dict(genes_of_scaffold_gpal.items() + genes_of_scaffold_gros.items())
	regions_of_scaffolds_d = dict(regions_of_scaffolds_gpal.items() + regions_of_scaffolds_gros.items())
	data = parseAnchorpoints(anchorpoints_f)
	parseBaseclusters(baseclusters_f)

	dict_of_connections = {}
	for idx, multiplicon in data.items():
		if not multiplicon.gros_scaffold in dict_of_connections:
			dict_of_connections[multiplicon.gros_scaffold] = {}
		dict_of_connections[multiplicon.gros_scaffold][multiplicon.gpal_scaffold] = multiplicon.counts['anchorpoints']
		if not multiplicon.gpal_scaffold in dict_of_connections:
			dict_of_connections[multiplicon.gpal_scaffold] = {}
		dict_of_connections[multiplicon.gpal_scaffold][multiplicon.gros_scaffold] = multiplicon.counts['anchorpoints']

	graph = nx.MultiGraph()
	graph.length = {}
	graph.genes = {}

	for scaffold, connections in dict_of_connections.items():
		graph.add_node(scaffold)
		graph.length[scaffold] = lengths[scaffold]
		graph.genes[scaffold] = genes_of_scaffold_gpal.get(scaffold, genes_of_scaffold_gros.get(scaffold, 0))
		for connection, weight in connections.items():
			graph.add_node(connection)
			graph.add_edge(scaffold, connection, weight=weight)
			graph.genes[connection] = genes_of_scaffold_gpal.get(connection, genes_of_scaffold_gros.get(connection, 0))

	graph.edges()
	#pos = nx.graphviz_layout(graph,prog="neato")
	#pos = nx.spring_layout(graph, iterations=20)
	pos = graphviz_layout(graph,prog="neato")
	cluster_d = {}
	with open("parse_iadhore_cluster.txt", 'w') as fh:
		fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('idx', 'count_gpal_scaffolds', 'span_gpal_scaffolds', 'count_gros_scaffolds', 'span_gros_scaffolds', 'gpal_scaffolds', 'gros_scaffolds'))
		for idx, subgraph in enumerate(nx.connected_component_subgraphs(graph)):
			cluster = Cluster(idx)
			with open("parse_iadhore_cluster" + str(idx) + ".txt", 'w') as fh_2:
				for node in subgraph.nodes():
					if node.startswith('GROS'):
						fh_2.write("#%s\t%s\n" % (node, "\t".join(regions_of_scaffolds_d[node]['genes'])))
						for key in regions_of_scaffolds_d[node]:
							if not key == 'genes':
								fh_2.write("%s\t%s\n" % (key, "\t".join(regions_of_scaffolds_d[node][key])))
						cluster.gros_lengths.append(lengths[node])
						cluster.gros_nodes.append(node)
					else:
						fh_2.write("#%s\t%s\n" % (node, "\t".join(regions_of_scaffolds_d[node]['genes'])))
						for key in regions_of_scaffolds_d[node]:
							if not key == 'genes':
								fh_2.write("%s\t%s\n" % (key, "\t".join(regions_of_scaffolds_d[node][key])))
						cluster.gpal_lengths.append(lengths[node])
						cluster.gpal_nodes.append(node)
				fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (cluster.idx, len(cluster.gpal_lengths), sum(cluster.gpal_lengths), len(cluster.gros_lengths), sum(cluster.gros_lengths), ",".join(cluster.gpal_nodes), ",".join(cluster.gros_nodes)))

	fig = plt.figure(1,figsize=(200,200))
	edge_labels=dict([((u,v,),d['weight']) for u,v,d in graph.edges(data=True)])
	prefix_colour_d = {'GROS' : '#A0CBE2', 'path' : '#E69F00' }
	node_colours=[prefix_colour_d[u[0:4]] for u in graph.nodes()]
	edge_bbox = dict(boxstyle="square,pad=0.3", lw=2, fc = 'white')
	nx.draw_networkx_edge_labels(graph, pos, font_size = 20, edge_labels = edge_labels, label_pos=0.5, bbox=edge_bbox)
	nx.draw_networkx(graph, pos, node_size=[graph.genes[v]*200 for v in graph],font_size = 20, with_labels=True, node_color=node_colours, alpha = 0.7)
	nx.draw_networkx_edges(graph, pos, alpha = 0.4)
	nx.write_gml(graph, "synteny_clusters.gml")
	plt.savefig("synteny_clusters.pdf",dpi=50, format='pdf')

