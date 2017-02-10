# **tree2neo**
A tool to import Tree files produced  by FastTree to a Neo4j Graph database.

This tool does the following:
    
   * Starts up a Neo4j docker container
   * Downloads the reference [database](https://zenodo.org/record/252101#.WIHfgvF95hH) in the `/data` directory of the container
   * Gets TREE files from a provided directory 
   * Maps and loads the TREE data in Graph database
    

## Usage

**Clone this repository:**

```
$ git clone https://github.com/SANBI-SA/tree2neo.git
$ cd tree2neo
```

   * **Standalone :computer: :**
   
        ```
        $ virtualenv envname
        $ source envname/bin/activate
        $ pip install -r requirements.txt
        $ pip install --editable .
        $ tree2neo --help
        $ tree2neo init -d data/tree data/db/data
        ```
   * **Using docker/docker-compose :whale: :**
      
      ```
      $ docker-compose up --build -d
      ```