data = vcf.Reader(open(vcf_file))

dic = {}
    for record in data:
        key = [field['GT'] for field in record.samples]
        key = ''.join(key)
        if key not in dic:
            dic[key] = 1
        elif key in dic:
            value = dic[key] + 1
            dic[key] = value

df = pd.DataFrame.from_dict(dic, orient='index', columns = ['SNPs'])
