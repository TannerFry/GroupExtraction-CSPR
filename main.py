import sys

class CSPRparser:
    def __init__(self, inputFileName):
        self.multiSum = 0  # multitargetting sum taken from the previous version of make_graphs
        self.multiCount = 0  # multitargetting count taken from the previous version of make_graphs
        self.genome = ""  # genome name
        self.misc = ""  # anything from the misc line
        self.repeats = {}  # dictionary of the number of repeats. See the read_repeats function for more info
        self.seeds = {}  # dictionary of which chromesomes are repeats. See the read_repeats function for more info
        self.dec_tup_data = {}
        self.chromesomesSelectedList = list()
        self.popData = {}
        self.base_array_64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789=/"
        # file path variable
        self.fileName = inputFileName
        self.endo_info = dict()
        self.karystatsList = []
        self.seed_id_seq_pair ={}

    def endo_import(self):
        f = open("CASPERinfo")
        while True:
            line = f.readline()
            if line.startswith("ENDONUCLEASES"):
                break
        while True:
            line = f.readline()
            if line.startswith("-"):
                break
            else:
                myinfo = line.split(";")
                self.endo_info[myinfo[0]] = myinfo[1:]  # first is PAM list, second is seed length, third is tot length

    def int2nt(self, num):
        if num == 0:
            return 'A'
        elif num == 1:
            return 'T'
        elif num == 2:
            return 'C'
        elif num == 3:
            return 'G'
        else:
            return 'N'

    def nt2int(self,nt):
        if nt == 'A':
            return 0
        elif nt == 'T':
            return 1
        elif nt == 'C':
            return 2
        elif nt == 'G':
            return 3
        else:
            return 0

    def decompress64(self, base64seq, slength=0, toseq=False):
        base10seq = int()
        if isinstance(base64seq, str):
            for i in range(len(base64seq)):
                power = len(base64seq) - (i+1)
                index = self.base_array_64.find(base64seq[i])
                if index != -1:
                    base10seq += index*pow(64, power)
        else:
            base10seq = base64seq
        if toseq:
            seq = str()
            number = base10seq
            while number >= 4:
                rem = number % 4
                number = int(number/4)
                seq += self.int2nt(rem)
            seq += self.int2nt(number)
            for i in range(len(seq), slength):
                seq += 'A'
            return seq
        else:
            return base10seq

    def decompress_csf_tuple(self, locseq, bool=False, endo="spCas9"):
        # Lookup endonuclease sequence lengths for parsing
        if(bool == False):
            mytuple = locseq[:-1].split(",")
        else:
            mytuple = locseq.split(",")
            front_seq = mytuple[3]

        loc = self.decompress64(mytuple[0])
        seq = mytuple[1]
        scr = self.decompress64(mytuple[2])
        strand = seq.find("+")
        if strand != -1:
            dira = "+"
            sequence = seq[:strand]
            pam = seq[strand+1:]
        else:
            seq = seq.split("-")
            sequence = seq[0]
            pam = seq[1]
            dira = "-"
        if bool:
            seqlength = int(self.endo_info[endo][2]) - int(
                self.endo_info[endo][1])  # gets the tail sequence length for processing repeats
        else:
            seqlength = int(self.endo_info[endo][2])  # gets the total sequence length
        pamlength = len(self.endo_info[endo][0].split(",")[0])  # gets the length of the primary PAM
        sequence = self.decompress64(sequence, seqlength, True)
        pam = self.decompress64(pam, pamlength, True)
        # The for loops fixes the problem of A's not being added to the end because they are removed on compression
        if(bool == True):
            sequence = sequence + front_seq
        return int(loc), str(sequence), pam, int(scr), dira, endo

    def read_first_lines(self):
        fileStream = open(self.fileName, 'r')
        self.genome = fileStream.readline()
        colonIndex = self.genome.find(':') + 2
        buffer1 = self.genome[colonIndex:]
        self.genome = buffer1
        buffer = fileStream.readline()
        self.misc = fileStream.readline()
        colonIndex = self.misc.find(':') + 2
        buffer1 = self.misc[colonIndex:]
        self.misc = buffer1
        colonIndex = buffer.find(':') + 2
        for i in range(colonIndex, len(buffer)):
            if buffer[i] == ',':
                bufferString1 = buffer[colonIndex:i]
                colonIndex = i + 1
                self.karystatsList.append(int(bufferString1))
        fileStream.close()

    def read_repeats(self, endoChoice):
        index = 0

        seedLength = int(self.endo_info[endoChoice][1])

        # clear what is already in there
        self.repeats.clear()
        self.seeds.clear()

        # only read the repeats section of the file
        fileStream = open(self.fileName, 'r')
        buf = fileStream.readline()
        while buf != "REPEATS\n":
            buf = fileStream.readline()
        split_info = fileStream.read().split('\n')
        fileStream.close()

        # parse the info now and store it in the correct dictionaries
        while (index + 1 < len(split_info)):
            seed = self.decompress64(split_info[index], slength=seedLength)
            repeat = split_info[index + 1].split("\t")
            self.repeats[seed] = 0
            self.seeds[seed] = []
            self.dec_tup_data[seed] = []
            for item in repeat:
                if item != "":
                    self.repeats[seed] += 1
                    sequence = item.split(',')
                    self.seeds[seed].append(sequence)
                    temp = sequence[1:4]
                    temp.append(str(self.decompress64(seed, toseq=True, slength=int(seedLength))))
                    string = ",".join(temp)
                    dec = self.decompress_csf_tuple(string, bool=True, endo=endoChoice)
                    dec = list(dec)
                    dec.append(sequence[0])
                    self.dec_tup_data[seed].append(dec)
                    self.multiSum += self.decompress64(sequence[3], slength=seedLength)
                    self.multiCount += 1
            index = index + 2

        i = 0
        for item in self.dec_tup_data:
            sorted_temp = []
            temp = []
            j = 0
            while j < len(self.dec_tup_data[item]):
                temp.append(self.dec_tup_data[item][j])
                while (True):
                    if j < len(self.dec_tup_data[item]) - 1:
                        if self.dec_tup_data[item][j][6] == self.dec_tup_data[item][j + 1][6]:
                            temp.append(self.dec_tup_data[item][j + 1])
                            j += 1
                        else:
                            temp.sort(key=lambda x: x[0])
                            for sub in temp:
                                sorted_temp.append(sub)
                            temp = []
                            break
                    else:
                        temp.sort(key=lambda x: x[0])
                        for sub in temp:
                            sorted_temp.append(sub)
                        temp = []
                        break
                j += 1
            self.dec_tup_data[item] = sorted_temp
            i += 1
            i += 1

    def chro_bar_data(self, input_seed, endo):
        dic_info = {}
        seqLength = int(self.endo_info[endo][1])
        for seed in self.seeds:
            temp = seed
            temp1 = str(self.decompress64(temp, slength=seqLength, toseq=True))
            self.seed_id_seq_pair[temp1] = seed
            dic_info[temp1] = {}
            for repeat in self.seeds[seed]:
                if repeat[0] in dic_info[temp1]:
                    dic_info[temp1][repeat[0]].append(self.decompress64(repeat[1]))
                else:
                    dic_info[temp1][repeat[0]] = [self.decompress64(repeat[1])]
        self.info=dic_info
        self.fill_Chromo_Text(input_seed, dic_info)

    def fill_Chromo_Text(self, seed, info, dist=0):
        self.chromo_length = self.karystatsList
        self.chromo_pos = {}
        self.seq_data = []
        self.export_data = []
        self.bar_coords = []
        for chromo in info[seed]:
            pos = []
            for position in info[seed][chromo]:
                self.seq_data.append(seed)
                test1 = position / self.chromo_length[int(chromo) - 1]
                test1 = int(test1 * 485)
                pos.append(test1)
            self.chromo_pos[chromo] = pos
        self.bar_coords.clear()  # clear bar_coords list before creating visual
        ind = 0
        r = 0
        i = 0
        index = 0
        for chromo in self.chromo_pos:
            self.chromo_pos[chromo] = sorted(self.chromo_pos[chromo])
            num = 0
            while num < len(self.chromo_pos[chromo]):
                found = False
                k = self.chromo_pos[chromo][num]
                temp = []
                temp.append(ind)
                temp.append(k + 40)
                temp.append((i * 25) + 3 + 10 * i)
                temp.append((i * 25) + 22 + 10 * i)
                ind_1 = ind
                seq_1 = str(self.seq_data[ind_1])
                seed_id_1 = self.seed_id_seq_pair[seq_1]
                temp_1 = self.dec_tup_data[seed_id_1]
                temp_1 = temp_1[ind]
                loc_1 = int(temp_1[0])
                for coords in self.bar_coords:
                    for sub_coords in coords:
                        ind_2 = sub_coords[0]
                        seq_2 = str(self.seq_data[ind_2])
                        seed_id_2 = self.seed_id_seq_pair[seq_2]
                        temp_2 = self.dec_tup_data[seed_id_2]
                        temp_2 = temp_2[ind_2]
                        loc_2 = int(temp_2[0])
                        if temp[1:4] == sub_coords[1:4] or (
                                loc_2 < (loc_1 + dist) and loc_2 > (loc_1 - dist) and temp[2:4] == sub_coords[2:4]):
                            found = True
                            exp = temp
                            exp.append(loc_1)
                            exp.append(seq_1)
                            exp.append(seed_id_1)
                            self.export_data[index - 1].append(exp)
                            break
                    if found == True:
                        break
                if found != True:
                    exp = temp
                    exp.append(loc_1)
                    exp.append(seq_1)
                    exp.append(seed_id_1)
                    self.export_data.append([exp])
                    index += 1
                ind += 1
                num += 1
            i = i + 1
            r += 1

    def export(self):
        groups = {}
        for coord in self.export_data:
            found = False
            for key in groups.keys():
                temp = []
                for obj in groups[key]:
                    for line in obj:
                        temp.append(line[4])
                ma = max(temp)
                mi = min(temp)
                if coord[0][1] == key or (coord[0][4] <= ma + 700 and coord[0][4] >= mi - 700):
                    groups[key].append(coord)
                    found = True
                    break
            if found == False:
                groups[coord[0][1]] = []
                groups[coord[0][1]].append(coord)
        for key in groups:
            print(key)
            print(groups[key])
            print('\n\n')


filename = str(sys.argv[1])
seed_file = str(sys.argv[2])

cspr_name = filename[filename.rfind('\\')+1:]
endo = cspr_name[cspr_name.find('_')+1:cspr_name.find('.')]


parser = CSPRparser(str(sys.argv[1]))
parser.endo_import()
parser.read_first_lines()
parser.read_repeats(endoChoice=endo)

f = open('export_chromo_data.csv','w+')
f.close()

seeds = []
with open(seed_file, 'r') as f:
    for line in f:
        line = line.strip('\n')
        seeds.append(int(line))

seqLength = int(parser.endo_info[endo][1])
seed = str(parser.decompress64(seeds[0], slength=seqLength, toseq=True))
print(seed)
parser.chro_bar_data(seed, endo)

for seed in seeds:
    seqLength = int(parser.endo_info[endo][1])
    dec_seed = str(parser.decompress64(seed, slength=seqLength, toseq=True))
    print(dec_seed)
    parser.chro_bar_data(dec_seed, endo)

    groups = {}
    for coord in parser.export_data:
        found = False
        for key in groups.keys():
            temp = []
            for obj in groups[key]:
                for line in obj:
                    temp.append(line[4])
            ma = max(temp)
            mi = min(temp)
            if coord[0][1] == key or (coord[0][4] <= ma + 700 and coord[0][4] >= mi - 700):
                groups[key].append(coord)
                found = True
                break
        if found == False:
            groups[coord[0][1]] = []
            groups[coord[0][1]].append(coord)

    with open('export_chromo_data.csv','a') as f:
        f.write(dec_seed+'\n')
        i = 1
        for key in groups.keys():
            seed_data = {}
            for obj in groups[key]:
                for line in obj:
                    seed = str(parser.decompress64(line[6], slength=seqLength, toseq=True))
                    if line[5] not in seed_data:
                        seed_data[line[5]] = {}
                    if (line[4],seed) not in seed_data[line[5]]:
                        seed_data[line[5]][(line[4],seed)] = 1
                    else:
                        seed_data[line[5]][(line[4],seed)] += 1
            for seq in seed_data:
                for loc in seed_data[seq]:
                    f.write('Group ' + str(i) + ',' + str(loc[0]) + ',' + str(seed_data[seq][loc]) + ',' + str(loc[1]) + '\n')
            i += 1
