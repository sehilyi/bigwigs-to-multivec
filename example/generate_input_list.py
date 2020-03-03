with open('./sample_input/input_list.txt', 'w') as url_file:
    for i in range(5):
        url_file.writelines("sample_output/" + str(i + 1) + "_cistrome_db.bw\n")