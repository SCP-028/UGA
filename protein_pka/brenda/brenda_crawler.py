import requests
from bs4 import BeautifulSoup

BASEURL = "https://www.brenda-enzymes.org"

session = requests.Session()
headers = {
    "Accept-Language": "en-US;q=0.9,en;q=0.8,zh-CN,zh;q=0.7,zh-TW;q=0.6",
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Connection": "keep-alive",
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/68.0.3440.106 Safari/537.36"
}
session.headers.update(headers)

url = f"{BASEURL}/all_enzymes.php"
r = session.get(url)
soup = BeautifulSoup(r.text)
hrefs = soup.find_all("a", string="long")
url = f"{BASEURL}/enzyme.php"
hrefs = [url + "?" + x["href"].split("?")[1] + "&organism[]=Homo+sapiens" for x in hrefs]

with open("brenda_urls.txt", "w") as f:
    f.write("\n".join(hrefs))
