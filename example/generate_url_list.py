import urllib.request, json 

# Download:
# wget -i input_list.txt

def request():
    
    """
    Get the list of BigWig file urls by looking into all Cistrome DB Browser pages.
    """

    default_url = "http://dc2.cistrome.org/api/main_filter_ng?allqc=true&cellinfos=all&completed=false&curated=false&factors=H3K27ac&keyword=h3k27ac&page={}&peakqc=false&run=false&species=Homo+sapiens"
    num_pages = -1

    with open('./example/input_list.txt', 'w') as url_list:

        # Get metadata (e.g., number of pages)
        with urllib.request.urlopen(default_url.format(1)) as url:
            json_result = json.loads(url.read().decode())
            num_pages = json_result["num_pages"]
        
        # Get lists
        for p in range(num_pages):
            print("Done processing page", p)
            with urllib.request.urlopen(default_url.format(p + 1)) as url:
                datasets = json_result["datasets"]

                for d in datasets:
                    cid = d["id"]
                    input_url = "http://dbtoolkit.cistrome.org/api_bigwig?sid={}".format(cid)
                    url_list.writelines(input_url + "\n")

if __name__ == "__main__":
	request()