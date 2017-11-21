configfile: "config.yaml"

rule all:
	input:
		['bg_seq/SE_R1-R5_{region_len}.seq'.format(region_len=x) for x in config['region_len']],
		expand('{AS_event}.{region_len}bp.k{kmer}.txt', AS_event=config['AS_event'], region_len=config['region_len'], kmer=config['kmer'])

rule get_bg_seq:
	input:
		genome = config['genome'],
		ASTA = config['ASTA']
	output:
		['bg_seq/SE_R1-R5_{region_len}.seq'.format(region_len=x) for x in config['region_len']]
	params:
		region_len = config['region_len']
	shell:
		'python script/get_R1-R5.py {input.genome} {input.ASTA} bg_seq {params.region_len}'

rule AS_kmer_enrich:
	input:
		bg_seq = 'bg_seq/SE_R1-R5_{region_len}.seq',
		events = '{AS_event}'
	output:
		'{AS_event}.{region_len}bp.k{kmer}.txt'
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/AS_kmer_enrich.R {input.bg_seq} {input.events} {output}'
