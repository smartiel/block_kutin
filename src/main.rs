extern crate clap;

use clap::{command, Arg};
use kutin_enumeration_rust::enumeration::generate_full_database;
use kutin_enumeration_rust::operators::Topology;
// use std::collections::HashMap;

// fn display_depth_distribution(database: &HashMap<u64, Vec<u64>>) {
//     let mut counts: HashMap<usize, usize> = HashMap::new();
//     let mut mdepth = 0;
//     for (key, hist) in database.iter() {
//         *counts.entry(hist.len()).or_insert(0) += 1;
//         if hist.len() > mdepth {
//             mdepth = hist.len();
//         }
//     }
//     for i in 0..(mdepth + 1) {
//         if counts.contains_key(&i) {
//             println!("{i}: {}", counts[&i]);
//         }
//     }
// }

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
    let topology = Topology::from_string(&topology_name);
    generate_full_database(
        &topology,
        &args.get_one("OUTPUT_FILE").unwrap(),
        (*args.get_one::<String>("STEP").unwrap())
            .parse::<i32>()
            .unwrap(),
    )
    .expect("Couldn't open the output file");
}
