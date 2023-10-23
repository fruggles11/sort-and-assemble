package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"net/http"
	"os"
	"path"
	"strings"
	"sync"
	"time"
)

func parseURLs(http_file string) ([]string, error) {

	var urls []string

	// Open the file
	file, err := os.Open(http_file)
	if err != nil {
		return urls, err
	}
	defer file.Close()

	// Create a new Scanner for the file
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		// Read each line and append it to the slice
		urls = append(urls, scanner.Text())
	}

	// Check for errors in the Scanner
	if err := scanner.Err(); err != nil {
		return urls, err
	}

	return urls, err
}

func downloadFile(url string) error {

	var resp *http.Response
	var err error

	// Retry up to 3 times
	for retries := 0; retries < 3; retries++ {
		resp, err = http.Get(url)
		if err == nil && resp.StatusCode == http.StatusOK {
			break // Success, exit the retry loop
		}
		if err != nil {
			fmt.Printf("Error fetching URL: %s. Retrying...\n", err)
		} else {
			fmt.Printf("Received unexpected status code: %d. Retrying...\n", resp.StatusCode)
		}
		time.Sleep(2 * time.Second) // Wait before retrying
	}

	// If the request failed, return the error
	if err != nil || resp.StatusCode != http.StatusOK {
		return err
	}

	// Parse the last part of the URL to get the file name
	fileName := path.Base(url)
	if fileName == "." || strings.Contains(fileName, "/") {
		return err
	}

	// Create a new file with the parsed name
	file, err := os.Create(fileName)
	if err != nil {
		return err
	}
	defer file.Close()

	// Copy the contents of the HTTP response to the file
	_, err = io.Copy(file, resp.Body)
	if err != nil {
		return err
	}

	return err
}

func main() {

	// Collect provided HTTP links from text file
	httpPtr := flag.String("http", "file-https.txt", "Text file with one file's HTTP address per line.")

	// Parse the command line arguments
	flag.Parse()

	// retrieve all the desired URLs
	urls, err := parseURLs(*httpPtr)
	if err != nil {
		log.Fatalf("Failed to scan and pull URLs from provided file: %s", err)
	}

	// Define the number of concurrent workers
	const numWorkers = 100

	// Create a channel to send UrlWithID structs to workers
	urlsChannel := make(chan string, numWorkers)

	// Create a wait group to wait for all workers to complete
	var wg sync.WaitGroup

	// Launch workers
	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for url := range urlsChannel {
				err := downloadFile(url)
				if err != nil {
					fmt.Printf("Error downloading URL %s: %v\n", url, err)
				}
			}
		}()
	}

	// Send urls structs to the url channel
	for _, url := range urls {
		urlsChannel <- url
	}
	close(urlsChannel)

	// Wait for all workers to complete
	wg.Wait()

	fmt.Println("Download completed.")

}
