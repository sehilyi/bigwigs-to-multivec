import urllib.request, json

# Download:
# wget -i input_urls.txt

def request():
    
    """
    Get the list of BigWig file urls by looking into all Cistrome DB Browser pages.
    """

    allqc = True
    factors = "H3K27ac" # or "all"
    keyword = "h3k27ac" # or ""
    default_url = "http://dc2.cistrome.org/api/main_filter_ng?allqc={}&cellinfos=all&completed=false&curated=false&factors={}&keyword={}&page={}&peakqc=false&run=false&species=Homo+sapiens"
    limit_num = 99999
    
    num_pages = -1
    url_list = []
    id_list = []

    # Get metadata (e.g., number of pages)
    with urllib.request.urlopen(default_url.format(allqc, factors, keyword, 1)) as url:
        json_result = json.loads(url.read().decode())
        num_pages = json_result["num_pages"]
    print("Total num of pages:", num_pages)
    
    # Get lists
    done_processing = False
    for p in range(num_pages):
        try:
            with urllib.request.urlopen(default_url.format(allqc, factors, keyword, p + 1)) as url:
                json_result = json.loads(url.read().decode())
                datasets = json_result["datasets"]

                for d in datasets:
                    cid = d["id"]
                    new_url = "http://dbtoolkit.cistrome.org/api_bigwig?sid={}".format(cid)
                    if new_url not in url_list:
                        url_list += [new_url]
                        id_list += [cid]
                        if len(url_list) >= limit_num:
                            done_processing = True
                    if done_processing:
                        break
        except:
            print("Error requesting", default_url.format(allqc, factors, keyword, p + 1))
        
        if(p % 100 is 1):
            print("Done processing page", p)

        if done_processing:
            print("Done processing.")
            break
    
    # Write input urls.
    with open('./sample_input/input_urls.txt', 'w') as f_out:
        # Write parameters.
        f_out.writelines("allqc={}\tfactors={}\tkeyword={}\tlimit_num={}\n".format(allqc, factors, keyword, limit_num))

        for u in url_list:
            f_out.writelines(u + "\n")
    
    # Write input urls.
    with open('./sample_input/input_list.txt', 'w') as f_out:
        for cid in id_list:
            f_out.writelines("api_bigwig?sid={}\n".format(cid))

if __name__ == "__main__":
	request()