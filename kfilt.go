package main

import (
	"bufio"
	"compress/gzip"
	"encoding/gob"
	"fmt"
	"io"
	"log"
	"os"
	"runtime"
	"strings"
	"sync"
	"sync/atomic"
	"time"

	"github.com/spf13/cobra"
	"github.com/cheggaaa/pb/v3"
)

const version = "1.0.0"

type BKNode struct {
	Kmer     string
	Children map[int]*BKNode
}

type BKTree struct {
	Root  *BKNode
	K     int
	Count int64
}

func NewBKTree(k int) *BKTree {
	return &BKTree{
		Root: &BKNode{Children: make(map[int]*BKNode)},
		K:    k,
	}
}

func (tree *BKTree) Insert(kmer string) {
	if tree.Root.Kmer == "" {
		tree.Root.Kmer = kmer
		tree.Count++
		return
	}
	if tree.Root.insert(kmer) {
		tree.Count++
	}
}

func (node *BKNode) insert(kmer string) bool {
	dist := hammingDistance(node.Kmer, kmer)
	if dist == 0 {
		return false
	}
	if child, exists := node.Children[dist]; exists {
		return child.insert(kmer)
	} else {
		node.Children[dist] = &BKNode{
			Kmer:     kmer,
			Children: make(map[int]*BKNode),
		}
		return true
	}
}

func (tree *BKTree) Search(kmer string, maxDist int) bool {
	if tree.Root.Kmer == "" {
		return false
	}
	return tree.Root.search(kmer, maxDist)
}

func (node *BKNode) search(kmer string, maxDist int) bool {
	dist := hammingDistance(node.Kmer, kmer)
	if dist <= maxDist {
		return true
	}
	minDist := dist - maxDist
	if minDist < 0 {
		minDist = 0
	}
	maxChildDist := dist + maxDist
	for childDist := minDist; childDist <= maxChildDist; childDist++ {
		if child, exists := node.Children[childDist]; exists {
			if child.search(kmer, maxDist) {
				return true
			}
		}
	}
	return false
}

func (tree *BKTree) SearchDetailed(kmer string, maxDist int) (bool, string, int) {
	if tree.Root.Kmer == "" {
		return false, "", -1
	}
	return tree.Root.searchDetailed(kmer, maxDist, "", 999999)
}

func (node *BKNode) searchDetailed(kmer string, maxDist int, bestMatch string, bestDist int) (bool, string, int) {
	dist := hammingDistance(node.Kmer, kmer)
	found := false
	if dist <= maxDist {
		found = true
		if dist < bestDist {
			bestDist = dist
			bestMatch = node.Kmer
		}
	}
	minDist := dist - maxDist
	if minDist < 0 {
		minDist = 0
	}
	maxChildDist := dist + maxDist
	for childDist := minDist; childDist <= maxChildDist; childDist++ {
		if child, exists := node.Children[childDist]; exists {
			childFound, childMatch, childDist := child.searchDetailed(kmer, maxDist, bestMatch, bestDist)
			if childFound {
				found = true
				if childDist < bestDist {
					bestDist = childDist
					bestMatch = childMatch
				}
			}
		}
	}
	return found, bestMatch, bestDist
}

func (tree *BKTree) Save(filename string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()
	gzWriter := gzip.NewWriter(file)
	defer gzWriter.Close()
	encoder := gob.NewEncoder(gzWriter)
	return encoder.Encode(tree)
}

func LoadBKTree(filename string) (*BKTree, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()
	gzReader, err := gzip.NewReader(file)
	if err != nil {
		return nil, err
	}
	defer gzReader.Close()
	var tree BKTree
	decoder := gob.NewDecoder(gzReader)
	if err := decoder.Decode(&tree); err != nil {
		return nil, err
	}
	return &tree, nil
}

var complementTable = [256]byte{
	'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
	'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
	'N': 'N', 'n': 'n',
}

func reverseComplement(seq string) string {
	rc := make([]byte, len(seq))
	for i := 0; i < len(seq); i++ {
		rc[len(seq)-1-i] = complementTable[seq[i]]
	}
	return string(rc)
}

func hammingDistance(s1, s2 string) int {
	if len(s1) != len(s2) {
		return len(s1)
	}
	dist := 0
	for i := 0; i < len(s1); i++ {
		if s1[i] != s2[i] {
			dist++
		}
	}
	return dist
}

func getReadBaseName(readName string) string {
	name := strings.TrimPrefix(strings.TrimPrefix(readName, "@"), ">")
	fields := strings.Fields(name)
	if len(fields) == 0 {
		return name
	}
	readID := fields[0]
	for _, suffix := range []string{"/1", "/2", ".1", ".2", "_1", "_2"} {
		if strings.HasSuffix(readID, suffix) {
			return readID[:len(readID)-2]
		}
	}
	if len(fields) >= 2 {
		lastField := fields[len(fields)-1]
		if strings.HasPrefix(lastField, "1:") || strings.HasPrefix(lastField, "2:") {
			return readID
		}
	}
	return readID
}

func areReadsPaired(read1Name, read2Name string) bool {
	base1 := getReadBaseName(read1Name)
	base2 := getReadBaseName(read2Name)
	return base1 == base2 && base1 != ""
}

func detectFormat(filename string) (string, error) {
	reader := openFile(filename)
	defer reader.Close()
	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)
	if scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "@") {
			return "fastq", nil
		} else if strings.HasPrefix(line, ">") {
			return "fasta", nil
		}
	}
	if err := scanner.Err(); err != nil {
		return "", err
	}
	return "", fmt.Errorf("could not detect file format (expected @ or > as first character)")
}

type ReadMatchInfo struct {
	ReadName        string
	IsR1            bool
	TotalKmers      int
	MatchedKmers    int
	BestHammingDist int
	Kept            bool
}

func writeVerbose(writer *bufio.Writer, info ReadMatchInfo) {
	readType := "SE"
	if info.IsR1 {
		readType = "R1"
	} else {
		readType = "R2"
	}
	kept := "FILTERED"
	if info.Kept {
		kept = "KEPT"
	}
	fmt.Fprintf(writer, "%s\t%s\t%d\t%d\t%d\t%s\n",
		strings.TrimPrefix(strings.TrimPrefix(info.ReadName, "@"), ">"),
		readType,
		info.TotalKmers,
		info.MatchedKmers,
		info.BestHammingDist,
		kept)
}

type FilterStats struct {
	TotalReads  int64
	KeptReads   int64
	PairedReads int64
	SingleReads int64
	StartTime   time.Time
	EndTime     time.Time
}

func buildCommand() *cobra.Command {
	var (
		kmerFile string
		output   string
		kmerSize int
	)
	cmd := &cobra.Command{
		Use:   "build",
		Short: "Build BK-tree index from k-mer list",
		Long:  "Build a BK-tree index from a k-mer list (e.g., meryl print output) for efficient approximate matching",
		RunE: func(cmd *cobra.Command, args []string) error {
			return runBuild(kmerFile, output, kmerSize)
		},
	}
	cmd.Flags().StringVarP(&kmerFile, "kmers", "k", "", "Input k-mer file (meryl print output)")
	cmd.Flags().StringVarP(&output, "output", "o", "bktree.idx", "Output BK-tree index file")
	cmd.Flags().IntVarP(&kmerSize, "kmer-size", "K", 31, "K-mer size")
	cmd.MarkFlagRequired("kmers")
	return cmd
}

func runBuild(kmerFile, output string, kmerSize int) error {
	log.Printf("Building BK-tree from %s (k=%d)...", kmerFile, kmerSize)
	totalLines, err := countLines(kmerFile)
	if err != nil {
		return fmt.Errorf("failed to count k-mers: %v", err)
	}
	tree := NewBKTree(kmerSize)
	file, err := os.Open(kmerFile)
	if err != nil {
		return err
	}
	defer file.Close()
	scanner := bufio.NewScanner(file)
	scanner.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)
	bar := pb.Full.Start64(totalLines)
	bar.Set(pb.Bytes, false)
	defer bar.Finish()
	count := int64(0)
	duplicates := int64(0)
	for scanner.Scan() {
		line := scanner.Text()
		if line == "" {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) > 0 {
			kmer := strings.ToUpper(fields[0])
			if len(kmer) == kmerSize {
				oldCount := tree.Count
				tree.Insert(kmer)
				if tree.Count == oldCount {
					duplicates++
				}
				count++
			}
		}
		bar.Increment()
	}
	if err := scanner.Err(); err != nil {
		return err
	}
	bar.Finish()
	log.Printf("Processed %d k-mers, %d unique, %d duplicates", count, tree.Count, duplicates)
	log.Printf("Saving index to %s...", output)
	if err := tree.Save(output); err != nil {
		return fmt.Errorf("failed to save index: %v", err)
	}
	info, _ := os.Stat(output)
	log.Printf("Index saved successfully (%.2f MB)", float64(info.Size())/1024/1024)
	return nil
}

func countLines(filename string) (int64, error) {
	file, err := os.Open(filename)
	if err != nil {
		return 0, err
	}
	defer file.Close()
	scanner := bufio.NewScanner(file)
	scanner.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)
	count := int64(0)
	for scanner.Scan() {
		if scanner.Text() != "" {
			count++
		}
	}
	return count, scanner.Err()
}

type Read struct {
	Name     string
	Sequence string
	Plus     string
	Quality  string
}

type ReadPair struct {
	R1          Read
	R2          Read
	IsPair      bool
	Matches     bool
	R1Matches   int
	R2Matches   int
	R1BestDist  int
	R2BestDist  int
}

type InputMode int

const (
	ModePairedSeparate InputMode = iota
	ModeSingleEnd
	ModeInterleaved
)

func filterCommand() *cobra.Command {
	var (
		input1        string
		input2        string
		interleaved   string
		index         string
		output        string
		outputFormat  string
		verboseFile   string
		minMatches    int
		hammingDist   int
		threads       int
		bothMatch     bool
		compress      bool
		compressLevel int
	)
	cmd := &cobra.Command{
		Use:   "filter",
		Short: "Filter reads based on k-mer matching",
		Long: `Filter FASTA/FASTQ reads that contain k-mers matching the BK-tree index.

Supports three input modes:
  1. Paired-end (separate files): -1 R1.fq -2 R2.fq
  2. Single-end: -1 reads.fq
  3. Interleaved (mixed paired/single): -I interleaved.fq

The tool automatically detects input format (FASTA or FASTQ).
Output format can be specified with --output-format (fasta or fastq).

For paired-end reads:
  - By default, keeps pair if combined matches >= threshold
  - With --both-match, requires BOTH reads to meet threshold individually

All formats support gzip compression with the -z flag.`,
		RunE: func(cmd *cobra.Command, args []string) error {
			return runFilter(input1, input2, interleaved, index, output, outputFormat,
				verboseFile, minMatches, hammingDist, threads, bothMatch, compress, compressLevel)
		},
	}
	cmd.Flags().StringVarP(&input1, "input1", "1", "", "Input FASTA/FASTQ R1 or single-end")
	cmd.Flags().StringVarP(&input2, "input2", "2", "", "Input FASTA/FASTQ R2 for paired-end")
	cmd.Flags().StringVarP(&interleaved, "interleaved", "I", "", "Interleaved FASTA/FASTQ (mixed paired/single)")
	cmd.Flags().StringVarP(&index, "index", "i", "bktree.idx", "BK-tree index file")
	cmd.Flags().StringVarP(&output, "output", "o", "", "Output file")
	cmd.Flags().StringVarP(&outputFormat, "output-format", "f", "fastq", "Output format: fasta or fastq")
	cmd.Flags().StringVarP(&verboseFile, "verbose", "v", "", "Verbose per-read output file (optional)")
	cmd.Flags().IntVarP(&minMatches, "min-matches", "n", 5, "Minimum matching k-mers required")
	cmd.Flags().IntVarP(&hammingDist, "hamming-dist", "m", 1, "Maximum Hamming distance")
	cmd.Flags().IntVarP(&threads, "threads", "t", runtime.NumCPU(), "Number of threads")
	cmd.Flags().BoolVarP(&bothMatch, "both-match", "B", false, "For pairs, require BOTH reads meet threshold")
	cmd.Flags().BoolVarP(&compress, "compress", "z", false, "Compress output with gzip")
	cmd.Flags().IntVarP(&compressLevel, "compress-level", "Z", 6, "Gzip compression level (1-9)")
	cmd.MarkFlagRequired("output")
	return cmd
}

func runFilter(input1, input2, interleaved, index, output, outputFormat, verboseFile string,
	minMatches, hammingDist, threads int, bothMatch, compress bool, compressLevel int) error {

	outputFormat = strings.ToLower(outputFormat)
	if outputFormat != "fastq" && outputFormat != "fasta" {
		return fmt.Errorf("invalid output format: %s (must be 'fastq' or 'fasta')", outputFormat)
	}

	var mode InputMode
	var inputFile1, inputFile2 string

	if interleaved != "" {
		if input1 != "" || input2 != "" {
			return fmt.Errorf("cannot specify both -I/--interleaved and -1/-2 flags")
		}
		mode = ModeInterleaved
		inputFile1 = interleaved
		log.Printf("Input mode: Interleaved (mixed paired/single-end)")
	} else if input1 != "" {
		if input2 != "" {
			mode = ModePairedSeparate
			inputFile1 = input1
			inputFile2 = input2
			log.Printf("Input mode: Paired-end (separate files)")
		} else {
			mode = ModeSingleEnd
			inputFile1 = input1
			log.Printf("Input mode: Single-end")
		}
	} else {
		return fmt.Errorf("must specify either -1/--input1 or -I/--interleaved")
	}

	format1, err := detectFormat(inputFile1)
	if err != nil {
		return fmt.Errorf("failed to detect format for %s: %v", inputFile1, err)
	}
	log.Printf("Detected input format: %s", strings.ToUpper(format1))
	log.Printf("Output format: %s", strings.ToUpper(outputFormat))

	var format2 string
	if mode == ModePairedSeparate {
		format2, err = detectFormat(inputFile2)
		if err != nil {
			return fmt.Errorf("failed to detect format for %s: %v", inputFile2, err)
		}
		if format1 != format2 {
			return fmt.Errorf("input files have different formats: %s vs %s", format1, format2)
		}
	}

	runtime.GOMAXPROCS(threads)

	log.Printf("Loading BK-tree from %s...", index)
	tree, err := LoadBKTree(index)
	if err != nil {
		return fmt.Errorf("failed to load index: %v", err)
	}
	log.Printf("Loaded index: k=%d, %d k-mers", tree.K, tree.Count)

	var outWriter *bufio.Writer
	outFile, err := os.Create(output)
	if err != nil {
		return err
	}
	defer outFile.Close()

	if compress {
		gzWriter, err := gzip.NewWriterLevel(outFile, compressLevel)
		if err != nil {
			return fmt.Errorf("failed to create gzip writer: %v", err)
		}
		defer gzWriter.Close()
		outWriter = bufio.NewWriter(gzWriter)
		log.Printf("Output compression: enabled (level %d)", compressLevel)
	} else {
		outWriter = bufio.NewWriter(outFile)
	}
	defer outWriter.Flush()

	var verboseWriter *bufio.Writer
	if verboseFile != "" {
		vFile, err := os.Create(verboseFile)
		if err != nil {
			return fmt.Errorf("failed to create verbose file: %v", err)
		}
		defer vFile.Close()
		verboseWriter = bufio.NewWriter(vFile)
		defer verboseWriter.Flush()
		verboseWriter.WriteString("read_id\tread_type\ttotal_kmers\tmatched_kmers\tbest_hamming_dist\tstatus\n")
		log.Printf("Verbose output: %s", verboseFile)
	}

	log.Println("Counting reads...")
	totalReads, err := countReadsByMode(inputFile1, inputFile2, mode, format1)
	if err != nil {
		log.Printf("Warning: couldn't count reads, progress bar will be indefinite: %v", err)
		totalReads = 0
	}

	if bothMatch && mode != ModeSingleEnd {
		log.Printf("Filtering mode: Both reads must meet threshold independently")
	} else {
		log.Printf("Filtering mode: Combined matches for pairs")
	}
	log.Printf("Filtering reads (min_matches=%d, hamming_dist=%d, threads=%d)...", 
		minMatches, hammingDist, threads)

	var stats *FilterStats
	switch mode {
	case ModePairedSeparate:
		r1Reader := openFile(inputFile1)
		defer r1Reader.Close()
		r2Reader := openFile(inputFile2)
		defer r2Reader.Close()
		stats = processReadsPaired(r1Reader, r2Reader, outWriter, verboseWriter, tree, 
			minMatches, hammingDist, threads, totalReads, bothMatch, format1, outputFormat)
	case ModeSingleEnd:
		reader := openFile(inputFile1)
		defer reader.Close()
		stats = processReadsSingle(reader, outWriter, verboseWriter, tree, 
			minMatches, hammingDist, threads, totalReads, format1, outputFormat)
	case ModeInterleaved:
		reader := openFile(inputFile1)
		defer reader.Close()
		stats = processReadsInterleavedMixed(reader, outWriter, verboseWriter, tree, 
			minMatches, hammingDist, threads, totalReads, bothMatch, format1, outputFormat)
	}

	log.Printf("\nFiltering complete!")
	log.Printf("Total reads: %d", stats.TotalReads)
	if mode == ModeInterleaved {
		log.Printf("  Paired reads: %d", stats.PairedReads)
		log.Printf("  Single reads: %d", stats.SingleReads)
	}
	log.Printf("Kept reads: %d (%.2f%%)", stats.KeptReads, 
		float64(stats.KeptReads)*100.0/float64(stats.TotalReads))
	log.Printf("Filtered reads: %d (%.2f%%)", stats.TotalReads-stats.KeptReads,
		float64(stats.TotalReads-stats.KeptReads)*100.0/float64(stats.TotalReads))
	log.Printf("Elapsed time: %.2fs", stats.EndTime.Sub(stats.StartTime).Seconds())
	log.Printf("Throughput: %.2f reads/sec", 
		float64(stats.TotalReads)/stats.EndTime.Sub(stats.StartTime).Seconds())

	return nil
}

func countReadsByMode(input1, input2 string, mode InputMode, format string) (int64, error) {
	switch mode {
	case ModePairedSeparate:
		count, err := countReads(input1, format)
		return count, err
	case ModeSingleEnd, ModeInterleaved:
		return countReads(input1, format)
	}
	return 0, nil
}

func openFile(filename string) io.ReadCloser {
	file, err := os.Open(filename)
	if err != nil {
		log.Fatal(err)
	}
	if strings.HasSuffix(filename, ".gz") {
		gzReader, err := gzip.NewReader(file)
		if err != nil {
			file.Close()
			log.Fatal(err)
		}
		return gzReader
	}
	return file
}

func countReads(filename string, format string) (int64, error) {
	reader := openFile(filename)
	defer reader.Close()
	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)
	count := int64(0)
	if format == "fastq" {
		for scanner.Scan() {
			line := scanner.Text()
			if strings.HasPrefix(line, "@") {
				count++
			}
		}
	} else {
		for scanner.Scan() {
			line := scanner.Text()
			if strings.HasPrefix(line, ">") {
				count++
			}
		}
	}
	return count, scanner.Err()
}

func parseFastaRecord(scanner *bufio.Scanner) (Read, bool) {
	var read Read
	if !scanner.Scan() {
		return read, false
	}
	read.Name = scanner.Text()
	if strings.HasPrefix(read.Name, ">") {
		read.Name = "@" + read.Name[1:]
	}
	if !scanner.Scan() {
		return read, false
	}
	read.Sequence = scanner.Text()
	read.Plus = "+"
	read.Quality = strings.Repeat("I", len(read.Sequence))
	return read, true
}

func parseFastqRecord(scanner *bufio.Scanner) (Read, bool) {
	var read Read
	if !scanner.Scan() {
		return read, false
	}
	read.Name = scanner.Text()
	if !scanner.Scan() {
		return read, false
	}
	read.Sequence = scanner.Text()
	if !scanner.Scan() {
		return read, false
	}
	read.Plus = scanner.Text()
	if !scanner.Scan() {
		return read, false
	}
	read.Quality = scanner.Text()
	return read, true
}

func parseRecord(scanner *bufio.Scanner, format string) (Read, bool) {
	if format == "fasta" {
		return parseFastaRecord(scanner)
	}
	return parseFastqRecord(scanner)
}

func writeReadFastq(w *bufio.Writer, read Read) {
	w.WriteString(read.Name)
	w.WriteByte('\n')
	w.WriteString(read.Sequence)
	w.WriteByte('\n')
	w.WriteString(read.Plus)
	w.WriteByte('\n')
	w.WriteString(read.Quality)
	w.WriteByte('\n')
}

func writeReadFasta(w *bufio.Writer, read Read) {
	name := read.Name
	if strings.HasPrefix(name, "@") {
		name = ">" + name[1:]
	}
	w.WriteString(name)
	w.WriteByte('\n')
	w.WriteString(read.Sequence)
	w.WriteByte('\n')
}

func writeRead(w *bufio.Writer, read Read, format string) {
	if format == "fasta" {
		writeReadFasta(w, read)
	} else {
		writeReadFastq(w, read)
	}
}

func countMatchingKmersDetailed(seq string, tree *BKTree, maxHamming int) (int, int, int) {
	k := tree.K
	if len(seq) < k {
		return 0, 0, 999
	}
	matches := 0
	totalKmers := 0
	bestDist := 999
	seqUpper := strings.ToUpper(seq)
	rcSeq := reverseComplement(seqUpper)
	for i := 0; i <= len(seqUpper)-k; i++ {
		kmer := seqUpper[i : i+k]
		totalKmers++
		found, _, dist := tree.SearchDetailed(kmer, maxHamming)
		if found {
			matches++
			if dist < bestDist {
				bestDist = dist
			}
		}
	}
	for i := 0; i <= len(rcSeq)-k; i++ {
		kmer := rcSeq[i : i+k]
		totalKmers++
		found, _, dist := tree.SearchDetailed(kmer, maxHamming)
		if found {
			matches++
			if dist < bestDist {
				bestDist = dist
			}
		}
	}
	if bestDist == 999 {
		bestDist = -1
	}
	return totalKmers, matches, bestDist
}

func processReadsPaired(r1Reader, r2Reader io.ReadCloser, outWriter, verboseWriter *bufio.Writer,
	tree *BKTree, minMatches, maxHamming int, threads int, totalReads int64, bothMatch bool, 
	inputFormat, outputFormat string) *FilterStats {

	stats := &FilterStats{StartTime: time.Now()}
	readChan := make(chan ReadPair, threads*4)
	resultChan := make(chan ReadPair, threads*4)

	var bar *pb.ProgressBar
	if totalReads > 0 {
		bar = pb.Full.Start64(totalReads)
		bar.Set(pb.Bytes, false)
		defer bar.Finish()
	}

	var writerWg sync.WaitGroup
	writerWg.Add(1)
	go func() {
		defer writerWg.Done()
		for pair := range resultChan {
			atomic.AddInt64(&stats.TotalReads, 1)
			atomic.AddInt64(&stats.PairedReads, 1)

			if verboseWriter != nil {
				totalKmers1 := 0
				if len(pair.R1.Sequence) >= tree.K {
					totalKmers1 = (len(pair.R1.Sequence) - tree.K + 1) * 2
				}
				totalKmers2 := 0
				if len(pair.R2.Sequence) >= tree.K {
					totalKmers2 = (len(pair.R2.Sequence) - tree.K + 1) * 2
				}
				writeVerbose(verboseWriter, ReadMatchInfo{
					ReadName:        pair.R1.Name,
					IsR1:            true,
					TotalKmers:      totalKmers1,
					MatchedKmers:    pair.R1Matches,
					BestHammingDist: pair.R1BestDist,
					Kept:            pair.Matches,
				})
				writeVerbose(verboseWriter, ReadMatchInfo{
					ReadName:        pair.R2.Name,
					IsR1:            false,
					TotalKmers:      totalKmers2,
					MatchedKmers:    pair.R2Matches,
					BestHammingDist: pair.R2BestDist,
					Kept:            pair.Matches,
				})
			}

			if pair.Matches {
				atomic.AddInt64(&stats.KeptReads, 1)
				writeRead(outWriter, pair.R1, outputFormat)
				writeRead(outWriter, pair.R2, outputFormat)
			}
			if bar != nil {
				bar.Increment()
			}
		}
	}()

	var workerWg sync.WaitGroup
	for i := 0; i < threads; i++ {
		workerWg.Add(1)
		go func() {
			defer workerWg.Done()
			for pair := range readChan {
				_, r1Matches, r1BestDist := countMatchingKmersDetailed(pair.R1.Sequence, tree, maxHamming)
				_, r2Matches, r2BestDist := countMatchingKmersDetailed(pair.R2.Sequence, tree, maxHamming)

				pair.R1Matches = r1Matches
				pair.R2Matches = r2Matches
				pair.R1BestDist = r1BestDist
				pair.R2BestDist = r2BestDist

				if bothMatch {
					pair.Matches = r1Matches >= minMatches && r2Matches >= minMatches
				} else {
					pair.Matches = (r1Matches + r2Matches) >= minMatches
				}
				resultChan <- pair
			}
		}()
	}

	go func() {
		scanner1 := bufio.NewScanner(r1Reader)
		scanner1.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)
		scanner2 := bufio.NewScanner(r2Reader)
		scanner2.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)
		for {
			r1, ok1 := parseRecord(scanner1, inputFormat)
			if !ok1 {
				break
			}
			r2, ok2 := parseRecord(scanner2, inputFormat)
			if !ok2 {
				log.Println("Warning: R2 ended before R1")
				break
			}
			readChan <- ReadPair{R1: r1, R2: r2, IsPair: true}
		}
		close(readChan)
	}()

	workerWg.Wait()
	close(resultChan)
	writerWg.Wait()
	if bar != nil {
		bar.Finish()
	}
	stats.EndTime = time.Now()
	return stats
}

func processReadsSingle(reader io.ReadCloser, outWriter, verboseWriter *bufio.Writer,
	tree *BKTree, minMatches, maxHamming int, threads int, totalReads int64, 
	inputFormat, outputFormat string) *FilterStats {

	stats := &FilterStats{StartTime: time.Now()}
	readChan := make(chan ReadPair, threads*4)
	resultChan := make(chan ReadPair, threads*4)

	var bar *pb.ProgressBar
	if totalReads > 0 {
		bar = pb.Full.Start64(totalReads)
		bar.Set(pb.Bytes, false)
		defer bar.Finish()
	}

	var writerWg sync.WaitGroup
	writerWg.Add(1)
	go func() {
		defer writerWg.Done()
		for pair := range resultChan {
			atomic.AddInt64(&stats.TotalReads, 1)
			atomic.AddInt64(&stats.SingleReads, 1)

			if verboseWriter != nil {
				totalKmers := 0
				if len(pair.R1.Sequence) >= tree.K {
					totalKmers = (len(pair.R1.Sequence) - tree.K + 1) * 2
				}
				writeVerbose(verboseWriter, ReadMatchInfo{
					ReadName:        pair.R1.Name,
					IsR1:            true,
					TotalKmers:      totalKmers,
					MatchedKmers:    pair.R1Matches,
					BestHammingDist: pair.R1BestDist,
					Kept:            pair.Matches,
				})
			}

			if pair.Matches {
				atomic.AddInt64(&stats.KeptReads, 1)
				writeRead(outWriter, pair.R1, outputFormat)
			}
			if bar != nil {
				bar.Increment()
			}
		}
	}()

	var workerWg sync.WaitGroup
	for i := 0; i < threads; i++ {
		workerWg.Add(1)
		go func() {
			defer workerWg.Done()
			for pair := range readChan {
				_, matches, bestDist := countMatchingKmersDetailed(pair.R1.Sequence, tree, maxHamming)
				pair.R1Matches = matches
				pair.R1BestDist = bestDist
				pair.Matches = matches >= minMatches
				resultChan <- pair
			}
		}()
	}

	go func() {
		scanner := bufio.NewScanner(reader)
		scanner.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)
		for {
			r1, ok := parseRecord(scanner, inputFormat)
			if !ok {
				break
			}
			readChan <- ReadPair{R1: r1, IsPair: false}
		}
		close(readChan)
	}()

	workerWg.Wait()
	close(resultChan)
	writerWg.Wait()
	if bar != nil {
		bar.Finish()
	}
	stats.EndTime = time.Now()
	return stats
}

func processReadsInterleavedMixed(reader io.ReadCloser, outWriter, verboseWriter *bufio.Writer,
	tree *BKTree, minMatches, maxHamming int, threads int, totalReads int64, bothMatch bool, 
	inputFormat, outputFormat string) *FilterStats {

	stats := &FilterStats{StartTime: time.Now()}
	readChan := make(chan ReadPair, threads*4)
	resultChan := make(chan ReadPair, threads*4)

	var bar *pb.ProgressBar
	if totalReads > 0 {
		bar = pb.Full.Start64(totalReads)
		bar.Set(pb.Bytes, false)
		defer bar.Finish()
	}

	var writerWg sync.WaitGroup
	writerWg.Add(1)
	go func() {
		defer writerWg.Done()
		for pair := range resultChan {
			atomic.AddInt64(&stats.TotalReads, 1)
			if pair.IsPair {
				atomic.AddInt64(&stats.PairedReads, 1)
			} else {
				atomic.AddInt64(&stats.SingleReads, 1)
			}

			if verboseWriter != nil {
				totalKmers1 := 0
				if len(pair.R1.Sequence) >= tree.K {
					totalKmers1 = (len(pair.R1.Sequence) - tree.K + 1) * 2
				}
				writeVerbose(verboseWriter, ReadMatchInfo{
					ReadName:        pair.R1.Name,
					IsR1:            true,
					TotalKmers:      totalKmers1,
					MatchedKmers:    pair.R1Matches,
					BestHammingDist: pair.R1BestDist,
					Kept:            pair.Matches,
				})

				if pair.IsPair {
					totalKmers2 := 0
					if len(pair.R2.Sequence) >= tree.K {
						totalKmers2 = (len(pair.R2.Sequence) - tree.K + 1) * 2
					}
					writeVerbose(verboseWriter, ReadMatchInfo{
						ReadName:        pair.R2.Name,
						IsR1:            false,
						TotalKmers:      totalKmers2,
						MatchedKmers:    pair.R2Matches,
						BestHammingDist: pair.R2BestDist,
						Kept:            pair.Matches,
					})
				}
			}

			if pair.Matches {
				atomic.AddInt64(&stats.KeptReads, 1)
				writeRead(outWriter, pair.R1, outputFormat)
				if pair.IsPair {
					writeRead(outWriter, pair.R2, outputFormat)
				}
			}
			if bar != nil {
				bar.Increment()
			}
		}
	}()

	var workerWg sync.WaitGroup
	for i := 0; i < threads; i++ {
		workerWg.Add(1)
		go func() {
			defer workerWg.Done()
			for pair := range readChan {
				_, r1Matches, r1BestDist := countMatchingKmersDetailed(pair.R1.Sequence, tree, maxHamming)
				pair.R1Matches = r1Matches
				pair.R1BestDist = r1BestDist

				if pair.IsPair {
					_, r2Matches, r2BestDist := countMatchingKmersDetailed(pair.R2.Sequence, tree, maxHamming)
					pair.R2Matches = r2Matches
					pair.R2BestDist = r2BestDist

					if bothMatch {
						pair.Matches = r1Matches >= minMatches && r2Matches >= minMatches
					} else {
						pair.Matches = (r1Matches + r2Matches) >= minMatches
					}
				} else {
					pair.Matches = r1Matches >= minMatches
				}
				resultChan <- pair
			}
		}()
	}

	go func() {
		scanner := bufio.NewScanner(reader)
		scanner.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)
		var previousRead *Read
		for {
			currentRead, ok := parseRecord(scanner, inputFormat)
			if !ok {
				if previousRead != nil {
					readChan <- ReadPair{R1: *previousRead, IsPair: false}
				}
				break
			}
			if previousRead == nil {
				previousRead = &currentRead
				continue
			}
			if areReadsPaired(previousRead.Name, currentRead.Name) {
				readChan <- ReadPair{R1: *previousRead, R2: currentRead, IsPair: true}
				previousRead = nil
			} else {
				readChan <- ReadPair{R1: *previousRead, IsPair: false}
				previousRead = &currentRead
			}
		}
		close(readChan)
	}()

	workerWg.Wait()
	close(resultChan)
	writerWg.Wait()
	if bar != nil {
		bar.Finish()
	}
	stats.EndTime = time.Now()
	return stats
}

func versionCommand() *cobra.Command {
	return &cobra.Command{
		Use:   "version",
		Short: "Print version information",
		Run: func(cmd *cobra.Command, args []string) {
			fmt.Printf("kfilt version %s\n", version)
			fmt.Printf("Go version: %s\n", runtime.Version())
			fmt.Printf("OS/Arch: %s/%s\n", runtime.GOOS, runtime.GOARCH)
		},
	}
}

func main() {
	rootCmd := &cobra.Command{
		Use:   "kfilt",
		Short: "Efficient k-mer based read filtering using BK-trees",
		Long: `kfilt: Fast and efficient k-mer based read filtering

This tool filters FASTA/FASTQ reads based on k-mer matching with support for
Hamming distance and reverse complement checking. It uses BK-trees for
efficient approximate k-mer matching.

Workflow:
  1. Build an index from k-mer list (e.g., from meryl difference)
  2. Filter multiple FASTA/FASTQ samples using the same index

The BK-tree index allows for efficient reuse across many samples without
rebuilding, making it ideal for processing large cohorts.`,
	}
	rootCmd.CompletionOptions.DisableDefaultCmd = true
	rootCmd.AddCommand(buildCommand())
	rootCmd.AddCommand(filterCommand())
	rootCmd.AddCommand(versionCommand())
	if err := rootCmd.Execute(); err != nil {
		os.Exit(1)
	}
}

