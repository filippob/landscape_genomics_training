#### support functions for the SMARTER-backend API

## 1) function to return the token needed for authenitcation
get_smarter_token <- function(base_url, username , password) {
  
  # username = readline(prompt = "Username ? "), password = getPass::getPass("Password ? ") # to use to prompt user/pwd request
  auth_url <-
    httr::modify_url(base_url, path = "/smarter-api/auth/login")
  
  resp <-
    POST(
      auth_url,
      body = list(username = username, password = password),
      encode = "json"
    )
  
  # this will read a JSON by default
  data <- httr::content(resp)
  
  # returning only the token as a string
  return(data$token)
}


## 2) function to read the URL
## (this function is used in get_smarter_data())
read_url <- function(url, token, query = list()) {
  # in this request, we add the token to the request header section
  resp <-
    GET(url, query = query, add_headers(Authorization = paste("Bearer", token)))
  
  # check errors: SMARTER-backend is supposed to return JSON objects
  if (http_type(resp) != "application/json") {
    stop("API did not return json", call. = FALSE)
  }
  
  # parse a JSON response. fromJSON to flatten results
  parsed <-
    jsonlite::fromJSON(
      content(resp, "text", encoding = "utf-8"),
      flatten = TRUE
    )
  
  # deal with API errors: not "200 Ok" status
  if (http_error(resp)) {
    stop(
      sprintf(
        "SMARTER API returned an error [%s]: '%s'",
        status_code(resp),
        parsed$message
      ),
      call. = FALSE
    )
  }
  
  return(parsed)
}

## 3) get SMARTER data
get_smarter_data <- function(url, base_url, token, query = list()) {
  # do the request and parse data with our function
  parsed <- read_url(url, token, query)
  
  # track results in df
  results <- parsed$items
  
  # check for pagination
  while (!is.null(parsed$`next`)) {
    # append next value to base url
    next_url <- httr::modify_url(base_url, path = parsed$`next`)
    
    # query arguments are already in url: get next page
    parsed <- read_url(next_url, token)
    
    # append new results to df. Deal with different columns
    results <- dplyr::bind_rows(results, parsed$items)
  }
  
  # return an S3 obj with the data we got
  structure(list(
    content = parsed,
    url = url,
    results = results
  ),
  class = "smarter_api")
}


#### 4) GET BREEDS
get_smarter_breeds <- function(bsurl, token, query = list()) {
  
  # setting the URL endpoint
  url <- httr::modify_url(bsurl, path = "/smarter-api/breeds")
  
  # reading our data
  data <- get_smarter_data(url, base_url = bsurl, token, query)
  
  # returning only the results dataframe
  data$results
}

#### GET DATASETS
get_smarter_datasets <- function(base_url, token, query = list()) {
  # setting the URL endpoint
  url <- httr::modify_url(base_url, path = "/smarter-api/datasets")
  
  # reading our data
  data <- get_smarter_data(url, base_url = base_url, token, query)
  
  # returning only the results dataframe
  data$results
}

#### GET SAMPLES
get_smarter_samples <- function(base_url, token, species, query = list()) {
  # mind that species is lowercase in endpoint url
  species <- tolower(species)
  
  url <- modify_url(base_url, path = sprintf("/smarter-api/samples/%s", species))
  
  data <- get_smarter_data(url, base_url = base_url, token, query)
  
  # returning only the results dataframe
  data$results
}



