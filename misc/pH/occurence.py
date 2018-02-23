import requests
import re
import time
from parsel import Selector
from urllib.parse import urljoin

headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.81 Safari/537.36',
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
    'Accept-Encoding': 'gzip, deflate, sdch, br',
    'Accept-Language': 'zh-CN,zh;q=0.8,en;q=0.6,zh-TW;q=0.4'
}
ss = requests.Session()
root_url = 'https://seer.cancer.gov/statfacts/more.html'
r = ss.get(root_url, headers=headers)
sel = Selector(r.text)
cancer_url = {}
for q in sel.xpath('//*[@id="content"]/div//ul/li'):
    url = q.xpath('./a/@href').extract_first()
    url = urljoin(root_url, url)
    cancer_name = q.xpath('./a/text()').extract_first()
    cancer_url[cancer_name] = url

cancer_num = {}
for cancer_name, url in cancer_url.items():
    try:
        r = ss.get(url, headers=headers)
        sel = Selector(r.text)
        q = sel.xpath('//*[@id="tog1"]/p[1]/text()').extract_first()
        case_num, death_num = re.findall("(\d+\.\d) per 100,000", q)
        cancer_num[cancer_name] = {
            'url': url,
            'case_num': case_num,
            'death_num': death_num
        }
        print('Finished cancer type: {}'.format(cancer_name))
        time.sleep(1)
    except Exception as e:
        print(e)

with open('cancer_occurence.csv', 'w') as f:
    f.write('cancer_name,case_num,death_num\n')
    for name, stat in cancer_num.items():
        print('{0},{1},{2}'.format(
            name, stat['case_num'], stat['death_num']), file=f)
