import argparse
import os
import subprocess
from helper_scripts.parseFastq import ParseFastQ
import pysam


class Mutation:
    
    def __init__(self, mutation, frequency, base_pos):
        self.mutation = mutation
        self.frequency = frequency
        self.base_pos = base_pos

    def __repr__(self):
        return f"{self.mutation}, {self.frequency} {self.base_pos}"

class ReportRow:

    def __init__(self, sample_name, color, total_reads, mutations: list[Mutation]):
        self.sample_name = sample_name
        self.color = color
        self.total_reads = total_reads
        self.mutations: list[Mutation] = mutations

    def __repr__(self):
        return f"name: {self.sample_name}, color: {self.color}, total reads: {self.total_reads}, mutations: {self.mutations}"

class ClinicalDataRow:

    def __init__(self, name, color, barcode):
        self.name = name
        self.color = color
        self.barcode = barcode

    def __repr__(self):
        return f"Name: {self.name} Barcode: {self.barcode} Color: {self.color}"


def main(fastq_file_name):
    process_and_trim(fastq_file_name)
    align()
    compress_alignments()
    report_list = pileup()
    create_report(report_list)


def create_report(report_list: list[ReportRow]):
    report_file = open("report.txt", "w")
    # pull over the report list containing all the mutations, their location, and names, and loop through them to make a report
    for report_row in report_list:
        print(f"Outputting report row {report_row.sample_name}")
        intro_str = f"Sample {report_row.sample_name} had a {report_row.color} mold, {report_row.total_reads} reads, and had "
        # append all the mutations
        mutations_str_list = []
        for mutation in report_row.mutations:
            percent_mutated = mutation.frequency / report_row.total_reads * 100
            mutations_str = f"{percent_mutated}% of the reads at position {mutation.base_pos} had the {mutation.mutation}"
            mutations_str_list.append(mutations_str)
        print("No mutations found")
        final_mutations_str = "no mutations"
        # if there are mutations append them
        if(len(mutations_str_list) > 0):
            final_mutations_str = ",".join(mutations_str_list)
        # output final report
        report_file.write(intro_str + final_mutations_str + ".\n")



def pileup():
    reference_file = open("dgorgon_reference.fa", "r")
    ref_seq = reference_file.readlines()[1].strip()
    report_list = []
    clinical_data = get_clinical_data()

    for clinical_data_row in clinical_data:
        samfile = pysam.AlignmentFile(f"sorted_bams/{clinical_data_row.name}.sorted.bam", "rb")
        #Since our reference only has a single sequence, we're going to pile up ALL of the reads. Usually you would do it in a specific region (such as chromosome 1, position 1023 to 1050 for example)
        for pileupcolumn in samfile.pileup():
            #use a dictionary to count up the bases at each position
            ntdict = {}
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    # get the base at the query postion
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    # get the reference base, and check if it equals the base at the query postion
                    ref_base = ref_seq[pileupread.query_position] if pileupread.query_position < len(ref_seq) else ""
                    if(base != ref_base):
                        # make a tuple with (frequency, query position) mapped to the base
                        if(base not in ntdict.keys()):
                            ntdict[base] = (1, pileupread.query_position)
                        else:
                            ntdict[base] = (ntdict[base][0], pileupread.query_position)
        # loop through the ntdict items looking for mutations, their frequency, and their base position 
        mutations: list[Mutation] = []
        for base, freq_base_pos in ntdict.items():
            mutation = Mutation(mutation=base, frequency=freq_base_pos[0], base_pos=freq_base_pos[1])
            mutations.append(mutation)
        print(f"Creating mutation report row for {clinical_data_row.name}")
        # create a report row that has the relevant clinical data, and the mutation data along with it.
        report_row = ReportRow(sample_name=clinical_data_row.name, color=clinical_data_row.color, total_reads=pileupcolumn.n, mutations=mutations)
        report_list.append(report_row)
        samfile.close()
    return report_list


def compress_alignments():
    print("Compressing and sort alignments into BAM files")
    clinical_data = get_clinical_data()

    os.makedirs("bams", exist_ok=True)
    os.makedirs("sorted_bams", exist_ok=True)

    # Create a bam file and sort it for each sam file
    for clinical_data_row in clinical_data:
        with open(f"bams/{clinical_data_row.name}.bam", "w") as output_file:
            try:
                # compress the sam file into a bam file
                compress_command = ["samtools", "view", "-bS", f"sams/{clinical_data_row.name}.sam"]
                subprocess.run(compress_command, check=True, stdout=output_file)
            except subprocess.CalledProcessError as e:
                print(f"Error executing command: {e}")
        # sort the bam file for better indexing
        sort_command = ["samtools", "sort", "-m", "100M", "-o", f"sorted_bams/{clinical_data_row.name}.sorted.bam", f"bams/{clinical_data_row.name}.bam"]
        subprocess.run(sort_command, check=True)
        index_command = ["samtools", "index", f"sorted_bams/{clinical_data_row.name}.sorted.bam"]
        subprocess.run(index_command, check=True)


def align():
    print("Aligning data...")
    clinical_data = get_clinical_data()
    # create the dgorgon_reference
    initial_align_command = ["bwa", "index", "dgorgon_reference.fa"]
    try:
        subprocess.run(initial_align_command, check=True)
        print("Initial align command executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        exit()

    os.makedirs("sams", exist_ok=True)

    # for each fastq file, align it to the reference
    for clinical_data_row in clinical_data:
        align_command = ["bwa", "mem", "dgorgon_reference.fa", f"fastqs/{clinical_data_row.name}_trimmed.fastq"]
        with open(f"sams/{clinical_data_row.name}.sam", "w") as output_file:
            try:
                subprocess.run(align_command, check=True, stdout=output_file)
                print("Align Command executed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error executing command: {e}")
                exit()


def process_and_trim(fastq_file_name):
    print("Processing and trimming data...")
    clinical_data = get_clinical_data()

    # for each clinical data row, process and trim it
    for clinical_data_row in clinical_data:
        print(f"Processing and trimming clinical data for {clinical_data_row.name}")
        fastq_file = ParseFastQ(fastq_file_name)

        # create the fastqs directory and trimmed files
        directory = 'fastqs'
        file_path = os.path.join(directory, f"{clinical_data_row.name}_trimmed.fastq")

        os.makedirs(directory, exist_ok=True)
        fastq_write_file = open(file_path, 'w')

        # for each fastq row, trim it and get sequence
        for fastq_row in fastq_file:
            if(fastq_row.seq_barcode == clinical_data_row.barcode):
                fastq_write_file.write(fastq_row.seq_header + "\n")
                fastq_write_file.write(fastq_row.seq_trimmed + "\n")
                fastq_write_file.write(fastq_row.qual_header + "\n")
                fastq_write_file.write(fastq_row.qual_trimmed + "\n")


# used to get clinical data in a structured format from the clinical data file
def get_clinical_data():
    clinical_data: list[ClinicalDataRow] = []
    try:
        clinical_data_file = open("harrington_clinical_data.txt", "r")
    except:
        print("ERROR: the clinical data file 'harrington_clinical_data.txt' must be in the same directory as this script")
        exit()

    for clinical_data_row in clinical_data_file.readlines()[1:]:
        clinical_data_row = clinical_data_row.split('\t')
        clinical_data.append(ClinicalDataRow(name=clinical_data_row[0].strip(), color=clinical_data_row[1].strip(), barcode=clinical_data_row[2].strip()))

    return clinical_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fastq", required=True, help="Place fastq inside here")
    args = parser.parse_args()
    main(args.fastq)
