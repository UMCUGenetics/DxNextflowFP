def flowcellLaneFromFastq(path) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab

    // parse first line of a FASTQ file (optionally gzip-compressed)
    // and return the flowcell id and lane number.
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven or eight fields)
    InputStream fileStream = new FileInputStream(path.toFile())
    InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)
    Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
    BufferedReader buffered = new BufferedReader(decoder)
    def line = buffered.readLine()
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(' ')[0].split(':')
    String machine
    int run_nr
    String fcid
    int lane

    machine = fields[0]
    run_nr = fields[1].toInteger()
    fcid = fields[2]
    lane = fields[3].toInteger()

    [fcid, lane, machine, run_nr]
}

def extractFastqPairFromDir(input, output) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    input = input.tokenize().collect{"$it/**_R1_*.fastq.gz"}
    analysis_id = output.split('/')[-1]
    Channel
    .fromPath(input, type:'file')
    .ifEmpty { error "No R1 fastq.gz files found in ${input}." }
    .filter { !(it =~ /.*Undetermined.*/) }
    .map { r1_path ->
        def fastq_files = [r1_path]
        def sample_id = r1_path.getSimpleName().split('_')[0]
        def r2_path = file(r1_path.toString().replace('_R1_', '_R2_'))
        if (r2_path.exists()) {
            fastq_files.add(r2_path)
        } else {
            exit 1, "R2 fastq.gz file not found: ${r2_path}."
        }
        def (flowcell, lane) = flowcellLaneFromFastq(r1_path)
        def rg_id = "${sample_id}_${flowcell}_${lane}"

        [['id': sample_id, 'rg_id': rg_id, 'flowcell': flowcell, 'analysis_id': analysis_id], fastq_files]
    }
}