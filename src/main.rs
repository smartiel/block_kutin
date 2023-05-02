extern crate clap;

use clap::{command, Arg};
use kutin_enumeration_rust::enumeration::{
    enumerate_step_1, enumerate_step_2, get_circuit_from_repr,
};
use kutin_enumeration_rust::operators::Topology;

fn main() {
    let args = command!()
        .version("0.0.1")
        .arg(
            Arg::new("TOPOLOGY")
                .help("The topology to run the enumeration on")
                .required(true),
        )
        .get_matches();
    let topology_name = args
        .get_one::<String>("TOPOLOGY")
        .expect("This should never fail :)");
    let topology = Topology::from_string(&topology_name);
    let db = enumerate_step_2(&topology);
    let operators = topology.get_all_operators();
    for repr in db.keys() {
        println!(
            "{}: {:?}",
            repr,
            get_circuit_from_repr(
                &db,
                &operators,
                repr,
                topology.nvertices,
                topology.nvertices
            )
        )
    }
}
