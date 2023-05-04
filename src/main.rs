extern crate block_kutin;
extern crate clap;

use block_kutin::enumeration::generate_database;
use clap::{command, Arg};

fn main() {
    let args = command!()
        .version("0.0.1")
        .arg(
            Arg::new("TOPOLOGY")
                .help("The topology to run the enumeration on")
                .required(true),
        )
        .arg(
            Arg::new("OUTPUT_FILE")
                .help("The file in which to dump the database")
                .required(true),
        )
        .arg(
            Arg::new("STEP")
                .help("The synthesis step (either 1 or 2, c.f. paper)")
                .required(true)
                .action(clap::ArgAction::Set),
        )
        .get_matches();
    let topology_name = args
        .get_one::<String>("TOPOLOGY")
        .expect("This should never fail :)");
    let file_name: &str = (*args.get_one::<&str>("OUTPUT_FILE").unwrap()).clone();
    generate_database(
        &topology_name,
        file_name,
        (*args.get_one::<String>("STEP").unwrap())
            .parse::<i32>()
            .unwrap(),
    )
    .expect("Couldn't open the output file");
}
