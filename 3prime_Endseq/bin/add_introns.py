#!/usr/local/bin/python

def add_introns_to_gtf(input_gtf_path, output_gtf_path):
    gene_transcripts = {}
    with open(input_gtf_path, 'r') as input_gtf:
        with open(output_gtf_path, 'w') as output_gtf:
            for line in input_gtf:
                if line.startswith('#'):
                    output_gtf.write(line)  # Write comments to output GTF as-is
                else:
                    fields = line.strip().split('\t')
                    feature_type = fields[2]
                    attributes = dict(item.strip().split(' ') for item in fields[8].strip().split(';') if item.strip())
                    transcript_id = attributes['transcript_id'].strip('"')
                    gene_id = attributes['gene_id'].strip('"')

                    key = (gene_id, transcript_id)
                    if key not in gene_transcripts:
                        gene_transcripts[key] = {'exons': [], 'other_features': []}

                    if feature_type == 'exon':
                        gene_transcripts[key]['exons'].append((fields, attributes))
                    else:
                        gene_transcripts[key]['other_features'].append(line)

            for (_, transcript_info) in gene_transcripts.items():
                exons = transcript_info['exons']
                exons.sort(key=lambda x: int(x[0][3]))  # Sort exons by start position

                # Write other features (CDS, UTRs) to output GTF
                for line in transcript_info['other_features']:
                    output_gtf.write(line)

                # Add introns between consecutive exons
                for i in range(1, len(exons)):
                    
                    (fields_exon, attributes_exon) = exons[i]
                    start_of_exon = fields_exon[3]
                   

                    (previous_exon_fields, _) = exons[i - 1]
                    end_of_previous_exon = previous_exon_fields[4]


                    intron_start = int(end_of_previous_exon) + 1
                    intron_end = int(start_of_exon) - 1
                    intron_number = i
                    fields[0] = fields_exon[0]
                    fields[1] = fields_exon[1]
                    fields[2] = 'intron'
                    fields[3] = str(intron_start)
                    fields[4] = str(intron_end)
                    fields[6] = fields_exon[6]
                    attributes['intron_number'] = str(intron_number)
                    

                    gene_id_current_exon = attributes_exon.get('gene_id', gene_id)
                    transcript_id_current_exon = attributes_exon.get('transcript_id', transcript_id)
                    attributes['uniq_trans_id'] = f"{transcript_id_current_exon[:-1]}_intron{intron_number}\""


                    # Write intron entry to output GTF
                    output_gtf.write('\t'.join(fields[:8]))
                    output_gtf.write(f"\tgene_id {gene_id_current_exon}; transcript_id {transcript_id_current_exon}; intron_number \"{intron_number}\"; uniq_trans_id {attributes['uniq_trans_id']};\n")

    print("Introns have been added to the GTF file.")

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Add introns to a GTF file with exons and other features.')
    parser.add_argument('input_gtf_file', help='Path to the input GTF file containing exons and other features.')
    parser.add_argument('output_gtf_file', help='Path to the output GTF file with added introns.')
    args = parser.parse_args()

    add_introns_to_gtf(args.input_gtf_file, args.output_gtf_file)
