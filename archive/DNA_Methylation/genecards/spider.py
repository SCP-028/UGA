import requests
import time

gene_list = ['AFAP1L2', 'ARSD', 'BAIAP2', 'BCL11A', 'CAST', 'CFLAR', 'CREM',
             'ECHDC2', 'GALNT6', 'GNAS', 'GRTP1', 'GSN', 'MLPH', 'MST1R',
             'NDRG2', 'PTK7', 'RGS14', 'RPS9', 'SEPT9', 'THNSL2']
url = 'http://www.genecards.org/GeneDownload/Enhancers'
params = {
    'GeneSymbol': ''
}
ss = requests.Session()
headers = {
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
    'Accept-Encoding': 'gzip, deflate, sdch',
    'Accept-Language': 'zh-CN,zh;q=0.8,en;q=0.6,zh-TW;q=0.4',
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.81 Safari/537.36'
}
for gene in gene_list:
    params['GeneSymbol'] = gene
    r = ss.get(url, headers=headers, params=params)
    with open('.\\{0}.txt'.format(gene), 'wb') as f:
        f.write(r.content)
    print('Finished querying gene: {0}'.format(gene))
    time.sleep(3)
