#!python
"""
Pull all ATP synthases.
"""
import requests
from parsel import Selector

df = []
r = requests.get('http://www.genenames.org/cgi-bin/genefamilies/set/644', headers= {"User-Agent": 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/60.0.3112.101 Safari/537.36',
         "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
         "Accept-Encoding": "gzip, deflate, sdch, br",
         "Accept-Language": "zh-CN,zh;q=0.8,en;q=0.6"})
sel = Selector(text=r.text)
for q in sel.xpath("//table//tr/td[1]/a/span"):
    df.append(q.xpath('text()').extract_first())
print('\n'.join(df))
