
#from collections import defaultdict

def decode_snpeff_annotations(vcf_file):
    annotation_details = []

    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  # Skip header lines
            columns = line.strip().split('\t')
            info_field = columns[7]
            info_parts = info_field.split(';')
            for part in info_parts:
                if part.startswith('ANN='):
                    annotations = part[4:].split(',')
                    for annotation in annotations:
                        details = annotation.split('|')
                        decoded_annotation = {
                            'allele': details[0],
                            'annotation': details[1],
                            'impact': details[2],
                            'gene_name': details[3],
                            'gene_id': details[4],
                            'feature_type': details[5],
                            'feature_id': details[6],
                            'transcript_biotype': details[7],
                            'rank_total': details[8],
                            'hgvs_c': details[9],
                            'hgvs_p': details[10],
                            'cdna_position_length': details[11],
                            'cds_position_length': details[12],
                            'aa_position_length': details[13],
                            'distance': details[14],
                            'error': details[15]
                        }
                        annotation_details.append(decoded_annotation)

    return annotation_details

def classify_annotations_by_impact(annotations):
    impact_categories = {
        'HIGH': [],
        'MODERATE': [],
        'LOW': [],
        'MODIFIER': []
    }

    for annotation in annotations:
        impact = annotation['impact']
        if impact in impact_categories:
            impact_categories[impact].append(annotation['annotation'])

    return impact_categories

vcf_file = '/Users/kenhsu/Dropbox/Mac/Desktop/311/高粱/master_thesis/WGSample/Out/30_vcf/qtlseq.ann.vcf'
annotations = decode_snpeff_annotations(vcf_file)
impact_categories = classify_annotations_by_impact(annotations)

# Print the categorized annotations
for impact, types in impact_categories.items():
    print(f"{impact} impact variants:")
    unique_types = set(types)
    for variant_type in unique_types:
        print(f"  - {variant_type} ({types.count(variant_type)} occurrences)")
    print()
