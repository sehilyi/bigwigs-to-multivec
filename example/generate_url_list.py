with open('./sample_input/url_list.txt', 'w') as url_file:
    for url in range(100):
        url_file.writelines("http://dbtoolkit.cistrome.org/api_bigwig?sid=" + str(url + 1) + "\n")

# Download:
# wget -i url_list.txt