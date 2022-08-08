import random, importlib, HTSeq, collections
import numpy as np

def set_bounds(txpt_id, db):
    
    CDS_bounds = []
    
    for feat in db.children(txpt_id, featuretype='CDS', order_by='start'):
        if len(CDS_bounds) == 0:
            CDS_bounds = [feat.start, feat.end]
        else:
            CDS_bounds = [ min([feat.start, CDS_bounds[0]]), max([feat.end, CDS_bounds[1]]) ]
                
    return CDS_bounds

def find_exons_for_genename(genename, db):
    
    five_regions = []
    cds_regions = []
    three_regions = []
    
    cds = set_bounds(genename, db)
    if len(cds) < 2:
        print(f"No CDS found for {genename}.")
        return {'5UTR': five_regions, 'CDS': cds_regions, '3UTR': three_regions}

    
    for feat in db.children(genename, featuretype='exon'):
        iv = ['', feat.start, feat.end]
        if feat.strand == '+' and iv[1] >= cds[1]: 
            three_regions.append(feat)
        elif feat.strand == '+' and iv[2] <= cds[0]:
            five_regions.append(feat)
        elif feat.strand == '-' and iv[1] <= cds[0]: 
            three_regions.append(feat)
        elif feat.strand == '-' and iv[2] >= cds[1]:
            five_regions.append(feat)
        elif (cds[0] <= iv[1] <= cds[1]) and (cds[0] <= iv[2] <= cds[1]):
            cds_regions.append(feat)
    
    return {'5UTR': five_regions, 'CDS': cds_regions, '3UTR': three_regions}



def plot_gene(iv, _ax, db, plot_this_txpt=None, xmin=0, xmax=None):

    #iv[0] = iv[0].split('chr')[-1]
    region = (str(iv[0]), iv[1], iv[2], str(iv[3]))

    wrote_gene = False
    txpt_to_tag = {}
    txpt_to_loc = {}
    gene_to_txpt = collections.defaultdict(set)
    print(f"Looking for gene in region:", region)
    for feat in db.region(region=region, featuretype=["exon"]):
        #print(f"Found {feat}")
        if 'tag' in feat.attributes:
            if 'basic' not in feat.attributes['tag']:
                #print(f"Skipped because 'tag' attribute was found without a 'basic' string.")
                continue
                
        if 'transcript_id' in feat.attributes:
            txpt = feat.attributes['transcript_id'][0]
            txpt_to_tag[txpt] = feat.attributes['tag']
            if txpt in txpt_to_loc:
                start, end = txpt_to_loc[txpt][0], txpt_to_loc[txpt][1]
                txpt_to_loc[txpt] = [min([feat.start, start]), max([feat.end, end])]
            else:
                txpt_to_loc[txpt] = [feat.start, feat.end]
            gene_to_txpt[feat.attributes['gene_name'][0]].add(txpt)
            
    print(txpt_to_loc)
    print(txpt_to_tag)
    print("gene_to_txpt:", gene_to_txpt)
    overlaps = {}
    for nameA, locA in txpt_to_loc.items():
        overlaps[nameA] = set([nameA])
        for nameB, locB in txpt_to_loc.items():
            if (locA[0] <= locB[0] <= locA[1]) or (locB[0] <= locA[0] <= locB[1]):
                overlaps[nameA].add(nameB)
    
    # Keep everything non-overlapping.
    txpt_to_plot = set([x for x in overlaps if len(overlaps[x])==0])
    
    print(f"Non-overlapping txpts to be plotted: {txpt_to_plot}")
    for gene, txpts in gene_to_txpt.items():
        has_tag = [txpt for txpt in txpts if (
            'appris_principal_1' in txpt_to_tag[txpt]) and (txpt not in txpt_to_plot)]

    if len(has_tag) == 1:
        txpt_to_plot.add(has_tag[0])
    elif len(has_tag) > 1:
        txpt_to_plot.add(sorted(has_tag, key=lambda x: float(x.split('ENST')[-1]))[0])
    else:
        txpt_to_plot.add(sorted(txpts, key=lambda x: float(x.split('ENST')[-1]))[0])
    
    if (plot_this_txpt is not None) and (plot_this_txpt in txpts):
        txpt_to_plot = set([plot_this_txpt])
        
    print(f"Drawing {txpt_to_plot}")
    for feat in db.region(region=region, featuretype=["exon"]):
        if feat.attributes['transcript_id'][0] not in txpt_to_plot:
            continue
        if 'gene_name' in feat.attributes:
            if 'SNOR' in feat.attributes['gene_name'][0] or "MIR" in feat.attributes['gene_name'][0]:
                continue
            #if wrote_gene and feat.attributes['gene_name'][0] != gene:
            #    continue
            if not wrote_gene:
                _ax.annotate(feat.attributes['gene_name'][0], xy=(0.05, 0.75), xycoords='axes fraction')
                wrote_gene = True
                gene = feat.attributes['gene_name'][0]
        left = max([iv[1], feat.start]) - iv[1]
        right = min([iv[2], feat.end]) - iv[1]
        #print(f"Drawing {left} to {right}.")
        lw = 3
        _ax.hlines(xmin=left, xmax=right, y=-1,  lw=lw, colors='k')

    for feat in db.region(region=(str(iv[0]), iv[1], iv[2], str(iv[3])), featuretype=["intron"]):
        if feat.attributes['transcript_id'][0] not in txpt_to_plot:
            continue
        left = max([iv[1], feat.start]) - iv[1]
        right = min([iv[2], feat.end]) - iv[1]
        lw = 1
        _ax.hlines(xmin=left, xmax=right, y=-1,  lw=lw, colors='k') 
        
    for feat in db.region(region=region, featuretype=['CDS']):
        if feat.attributes['transcript_id'][0] not in txpt_to_plot:
            continue
        left = max([iv[1], feat.start]) - iv[1]
        right = min([iv[2], feat.end]) - iv[1]
        lw = 6
        _ax.hlines(xmin=left, xmax=right, y=-1,  lw=lw, colors='k')

    #_ax.set_xlim((0, iv[2]-iv[1]))
    _ax.set_xlim((0, xmax))
    #_ax.axis('off')

    if wrote_gene:
        return gene
    return ''
