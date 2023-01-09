import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
from tkinter import StringVar
import _mypath
import InSiTE_module

# create the main window
root = tk.Tk()
root.title("Variable Input Form")

# create the label and text entry for the local filename
ttk.Label(root, text="Input File:").grid(row=0, column=0)
filename_entry = ttk.Entry(root)
filename_entry.grid(row=0, column=1)

# create the browse button for the local filename
def browse_filename():
    filepath = filedialog.askopenfilename()
    filename_entry.delete(0, tk.END)
    filename_entry.insert(0, filepath)

browse_button = ttk.Button(root, text="Browse", command=browse_filename)
browse_button.grid(row=0, column=2)

# create the dropdown menu for file type
ttk.Label(root, text="File type:").grid(row=1, column=0)
filetype_var = tk.StringVar(root)
filetype_var.set("fasta") # default value
filetype_menu = ttk.OptionMenu(root, filetype_var, "fasta", "fasta", "sam", "csv", "fastq")
filetype_menu.grid(row=1, column=1)

# create a function to update the dropdown menu to the selected value
def update_filetype(event):
    filetype_var.set(event)

# bind the update function to the dropdown menu
filetype_menu.bind("<<ComboboxSelected>>", update_filetype)


# set the default value of the string variable
string_var = StringVar()
string_var.set("CAMAGGTTGAAGAACACTG")

# create the label and text entry for the string input
string_label = ttk.Label(root, text="TPN Recognition Sequence")
string_label.grid(row=2, column=0)
string_entry = ttk.Entry(root, textvariable=string_var)
string_entry.grid(row=2, column=1, columnspan=2, sticky="W")

# create the tickbox for the remote option. help="get sequences from entrez server instead of 'LOCAL' TwoBit genome")
remote_var = tk.IntVar()
remote_checkbox = ttk.Checkbutton(root, text="Remote", variable=remote_var)
remote_checkbox.grid(row=2, column=3)

# create the label and text entry for the paired reads file
paired_reads_label = ttk.Label(root, text="Paired reads:")
paired_reads_label.grid(row=0, column=3)
paired_reads_entry = ttk.Entry(root)
paired_reads_entry.grid(row=0, column=4)

# create the browse button for the paired reads file
def browse_paired_reads():
    filepath = filedialog.askopenfilename()
    paired_reads_entry.delete(0, tk.END)
    paired_reads_entry.insert(0, filepath)

paired_reads_browse_button = ttk.Button(root, text="Browse", command=browse_paired_reads)
paired_reads_browse_button.grid(row=0, column=5)

# create a function to update the visibility of the paired reads field
def update_paired_reads_visibility(*args):
    filetype = filetype_var.get()
    if filetype == "fasta" or filetype == "fastq":
        paired_reads_label.grid()
        paired_reads_entry.grid()
        paired_reads_browse_button.grid()
    else:
        paired_reads_label.grid_remove()
        paired_reads_entry.grid_remove()
        paired_reads_browse_button.grid_remove()

# bind the update function to the file type dropdown menu
filetype_var.trace("w", update_paired_reads_visibility)

# initialize the visibility of the paired reads field
update_paired_reads_visibility()

# create the label and text entry for the vector. Help: FASTA file containing sequences te exclude from mapping,
# for example plasmids used in the experiment
vector = ttk.Label(root, text="Vector sequence")
vector.grid(row=3, column=0)
vector_entry = ttk.Entry(root)
vector_entry.grid(row=3, column=1)

# create the browse button for the vector filename
def browse_vector():
    filepath = filedialog.askopenfilename()
    vector_entry.delete(0, tk.END)
    vector_entry.insert(0, filepath)

browse_button = ttk.Button(root, text="Browse", command=browse_vector)
browse_button.grid(row=3, column=2)
# ap.add_argument('-r', '--rand_is', default=False, action='store_true',
#                 help='use random integration sites matched to given query set instead of actual query set')
# ap.add_argument('-n', '--rand_nt', default=False, action='store_true',
#                 help='use random sequences for mapping')

# ap.add_argument('--no_seqs', default=False, action='store_true',
#                 help='do not get sequences from either entrez or local TwoBit genome around locations indicated '
#                      'by genome_location_csv')

# ap.add_argument('--no_annotate', default=False, action='store_true',
#                 help='do not map insertion sites to genome annotations')
# ap.add_argument('--barcode', default='', metavar='NNNNN', help=' barcode sequence to trim off of reads')
# ap.add_argument('--samwindow', default=0, type=int, help='depreciated')  # 60

# #  ap.add_argument('--RC', default = None, metavar='file1 file2', help='XXX files to find & remove collisions, then reprocess')

# ap.add_argument('--primer3', default='', metavar='NNNNN', help="3' primer sequence to remove from reads")
# ap.add_argument('--trim5', default=0, type=int, help="additional (non-genomic) nts to trim off of 3' end of reads")
# ap.add_argument('--trim3', default=0, type=int,
#                 help="additional (non-genomic) nts to trim off of 5' end of reads")  # (starts with TA)")
# ap.add_argument('--feature', action='append', metavar='intron/exon/transcript/TSS/etc',
#                 help='feature names found in feature files to map reads to, e.g. "exon"')  # ['intron', 'exon',
# # 'codingexon','transcript','TSS'] #feature names found in feature files to map reads to
# ap.add_argument('--dist', action='append', metavar='True/False',
#                 help='weather to map distance of each read, or only whether reads overlap with feature. '
#                      'Same number of distance variables must be given as features.')
# ap.add_argument('--close', action='append', type=int, help='distance in bp to be considered close to feature')
# ap.add_argument('--chromosome_ids', metavar='/path/to/chromosomes.csv',
#                 default='./reference_datasets/chromosomes.csv')
# ap.add_argument('--bowtielocation', metavar='/path/to/bowtie2', default='/btapps/miniconda3/bin/bowtie2')
# ap.add_argument('--bowtieindex', metavar='/path/to/bowtieindex',
#                 default='./reference_datasets/genomes/GRCh38.fna.bowtie_index/'
#                         'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index')
# ap.add_argument('--weblogolocation', metavar='/path/to/weblogo', default='weblogo')  # in PATH
# ap.add_argument('--twobitlocation', metavar='/path/to/TwoBitToFa', default='./scripts/TwoBitToFa')
# ap.add_argument('--twobitgenomelocation', metavar='/path/to/genome.2bit',
#                 default='./reference_datasets/genomes/GRCh38.2bit')
# ap.add_argument('--annotations', action='append', metavar='/path/to/annotation_file.bed/gff/gtf',
#                 help='location of annotation file(s) (bed/gff/gtf), must be same number of files as features '
#                      'specified')  # ['./reference_datasets/annotations/gencode.v32.introns.bed',
# # './reference_datasets/annotations/gencode.v32.exons.bed',
# # './reference_datasets/annotations/gencode.v32.codingexons.bed',
# # './reference_datasets/annotations/gencode.v32.transcripts.gtf',
# # './reference_datasets/annotations/gencode.v32.transcripts.gtf']
# ap.add_argument('--supress_csv', default=False, action='store_true', help='do not output csv file')
# ap.add_argument('--supress_fasta', default=False, action='store_true', help='do not ouptput fasta file')
# ap.add_argument('--supress_logo', default=False, action='store_true', help='do not output logo')
# ap.add_argument('--append_summary', default=False, action='store_true',
#                 help='add summary metrics to summary.csv file')
# ap.add_argument('-v', '--verbose', default=False, action='store_true', help='verbose output and logging')
# args = ap.parse_args()


# create the hidden options frame
advanced_frame = tk.Frame(root)
advanced_frame.grid(row=4, column=0, columnspan=3)
advanced_frame.grid_remove()

# create the tickbox for the verbose option
verbose_var = tk.IntVar()
verbose_checkbox = ttk.Checkbutton(advanced_frame, text="Verbose", variable=verbose_var)
verbose_checkbox.grid(row=0, column=0)

# create the tickbox for the use random sequences option
use_random_sequences_var = tk.IntVar()
use_random_sequences_checkbox = ttk.Checkbutton(advanced_frame, text="Use random sequences", variable=use_random_sequences_var)
use_random_sequences_checkbox.grid(row=1, column=0)

# create the tickbox for the compress duplicate reads option
compress_duplicate_reads_var = tk.IntVar()
compress_duplicate_reads_var.set(1) # set to 1 (ticked) by default
compress_duplicate_reads_checkbox = ttk.Checkbutton(advanced_frame, text="Classify duplicate reads (+/- 1nt) as the same integration event", variable=compress_duplicate_reads_var)
compress_duplicate_reads_checkbox.grid(row=2, column=0)

# ap.add_argument('--lwindow', default=50, help='numebr of nucleotides upstream of integration site to return',
#                 type=int)  # window on either side of indicated nt location to return
# ap.add_argument('--rwindow', default=50, help='number of nucleotides downstream of integration site to return',
#                 type=int)  # window on either side of indicated nt location to return
# ap.add_argument('--min', default=25, type=int,
#                 help='minimum length of a read to try mapping. default (25) will usually avoid any false positives '
#                      'in read sets of 200k reads')


# create the advanced button to toggle the visibility of the hidden options frame
advanced_button = ttk.Button(root, text="Advanced")
advanced_button.grid(row=5, column=1)

def toggle_advanced_options():
    if advanced_frame.grid_info():
        advanced_frame.grid_remove()
    else:
        advanced_frame.grid()


advanced_button.config(command=toggle_advanced_options)

def submit():
    # get the values of the variables
    filename = filename_entry.get()
    filetype = filetype_var.get()
    string_input = string_entry.get()
    remote = remote_var.get()
    verbose = verbose_var.get()
    use_random_sequences = use_random_sequences_var.get()
    compress_duplicate_reads = compress_duplicate_reads_var.get()

    # close the main window
    root.destroy()

    # create a new window to display the results
    results_window = tk.Tk()
    results_window.title("Results")

    # display the results in the new window
    ttk.Label(results_window, text=f"Filename: {filename}").pack()
    ttk.Label(results_window, text=f"Filetype: {filetype}").pack()
    ttk.Label(results_window, text=f"String input: {string_input}").pack()
    ttk.Label(results_window, text=f"Remote: {remote}").pack()
    ttk.Label(results_window, text=f"Verbose: {verbose}").pack()
    ttk.Label(results_window, text=f"Use random sequences: {use_random_sequences}").pack()
    ttk.Label(results_window, text=f"Compress duplicate reads: {compress_duplicate_reads}").pack()

    # run the main loop for the results window
    results_window.mainloop()


# create the submit button and place it in the main window
submit_button = ttk.Button(root, text="Submit", command=submit)
submit_button.grid(row=6, column=1)

# run the main loop
root.mainloop()
