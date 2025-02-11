import argparse
#Example use is 
# python3 parseFastq.py --fastq /home/rbif/week5/hawkins_pooled_sequences.fastq

class FastQRow:

    def __init__(self, seq_header, seq_str, qual_header, qual_str):
        self.seq_header = seq_header
        self.seq_str = seq_str
        self.qual_header = qual_header
        self.qual_str = qual_str

    @property
    def seq_barcode(self):
        return self.seq_str[:5]
    
    @property
    def seq_trimmed(self):
        return self.trim_seq()[0]
    
    @property
    def qual_trimmed(self):
        return self.trim_seq()[1]

    def trim_seq(self):
        qual_beginning_trimmed = self.qual_str[5:]
        seq_beginning_trimmed = self.seq_str[5:]
        seq_trimmed = ""
        qual_trimmed = ""
        for i in range(2, len(qual_beginning_trimmed)):
            if(qual_beginning_trimmed[i - 2:i] in ["FD", "DF", "DD", "FF"]):
                return seq_trimmed, qual_trimmed
            qual_trimmed += qual_beginning_trimmed[i - 2]
            seq_trimmed += seq_beginning_trimmed[i - 2]
        return seq_trimmed, qual_trimmed

    def __repr__(self):
        return f"Seq Header: {self.seq_header}, Seq Barcode: {self.seq_barcode}"



################################################
# You can use this code and put it in your own script
class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self,filePath,headerSymbols=['@','+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...
 
        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'r')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
         
    def __iter__(self):
        return self
     
    def __next__(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)
         
        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 
         
        # ++++ Return fatsQ data as a FastQRow ++++
        return FastQRow(seq_header=elemList[0], seq_str=elemList[1], qual_header=elemList[2], qual_str=elemList[3])
##########################################################################



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fastq", required=True, help="Place fastq inside here")
    args = parser.parse_args()
    fastqfile = ParseFastQ(args.fastq)

    #A fastq read contains 4 lines
    for i, fastq_obj in enumerate(fastqfile):
        if(i == 1):
            break
        print(fastq_obj.seq_str)
        print(fastq_obj.qual_str)
        print('*'*10 + '==='*10 + '*' * 10)
        print(fastq_obj.seq_trimmed)
        print(fastq_obj.qual_trimmed)
        
        #Just an indicator showing the fastq "blocks"
        print('*'*10 + '==='*10 + '*' * 10)
