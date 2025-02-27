#BASE_DIR = "C:/Users/dshah35/OneDrive - St. Jude Children's Research Hospital/Desktop/files/Github/software/talent_tools_master"
#BASE_DIR = "."
import os
BASE_DIR = os.path.abspath(os.path.dirname(__file__))
GENOME_DIR = BASE_DIR + "/genome_data"
GENOME_FILE = GENOME_DIR + "/%s.fasta"
PROMOTEROME_DIR = BASE_DIR + "/promoterome_data"
PROMOTEROME_FILE = PROMOTEROME_DIR + "/%s.fasta"
RVD_SEQ_REGEX = r'^(?:(?:[ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY\*])[ _]+){11,30}(?:[ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY\*])$'
#DRUPAL_CALLBACK_URL = "http://talent.local/talent/jobcomplete/"
DRUPAL_CALLBACK_URL = "https://tale-nt.cac.cornell.edu/talent/jobcomplete/"
VALID_GENOME_ORGANISMS = ['drosophila_melanogaster', 'arabidopsis_thaliana', 'mus_musculus', 'oryza_sativa', 'caenorhabditis_elegans', 'danio_rerio', 'homo_sapiens', 'rattus_norvegicus', 'brachypodium_distachyon', 'solanum_lycopersicum', 'gasterosteus_aculeatus']
VALID_PROMOTEROME_ORGANISMS = ['drosophila_melanogaster', 'arabidopsis_thaliana', 'mus_musculus', 'oryza_sativa', 'oryza_sativa_msu7', 'oryza_sativa_msu6', 'oryza_sativa_rapdb_only', 'caenorhabditis_elegans', 'danio_rerio', 'homo_sapiens', 'rattus_norvegicus', 'brachypodium_distachyon', 'gasterosteus_aculeatus']
REDIS_SERVER_HOSTNAME = "localhost"
REDIS_SERVER_PORT = 6379
OFFTARGET_COUNTING_SIZE_LIMIT = 316000000