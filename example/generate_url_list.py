import urllib.request, json 
from urllib.error import HTTPError

# Download:
# wget -i input_list.txt

def request():
    
    """
    Get the list of BigWig file urls by looking into all Cistrome DB Browser pages.
    """

    default_url = "http://dc2.cistrome.org/api/main_filter_ng?allqc=false&cellinfos=all&completed=false&curated=false&factors=all&keyword=&page={}&peakqc=false&run=false&species=Homo+sapiens"
    num_pages = -1
    url_list = []

    # Get metadata (e.g., number of pages)
    with urllib.request.urlopen(default_url.format(1)) as url:
        json_result = json.loads(url.read().decode())
        num_pages = json_result["num_pages"]
    print("Total num of pages:", num_pages)
    
    # Get lists
    for p in range(num_pages):
        try:
            with urllib.request.urlopen(default_url.format(p + 1)) as url:
                datasets = json_result["datasets"]

                for d in datasets:
                    cid = d["id"]
                    new_url = "http://dbtoolkit.cistrome.org/api_bigwig?sid={}".format(cid)
                    if new_url not in url_list:
                        url_list += [new_url]
        except HTTPError as e:
            print("Error requesting", default_url.format(p + 1), e.read())
        if(p % 100 is 0):
            print("Done processing page", p)
                
    with open('./example/input_list.txt', 'w') as f_out:
        for u in url_list:
            f_out.writelines(u + "\n")

if __name__ == "__main__":
	request()