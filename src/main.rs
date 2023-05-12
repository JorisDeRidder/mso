
use rpassword::prompt_password;
use postgres::{Client, NoTls};




#[derive(Debug)]
struct RunInfo {
    runid: i32,
    runname: String
}




fn main() {

    // For security reasons, don't hardcode the user ID and password

    let username = prompt_password("Enter server username: ").expect("Failed to read username.");
    let password = prompt_password("Enter server password: ").expect("Failed to read password.");

    // Connect to the server

    let connection_string = format!("postgresql://{}:{}@gaiadb08i:55431/surveys", username, password);
    let mut client = Client::connect(&connection_string, NoTls).expect("Failed to connect to the server.");

    let query = "
        select 
            runid, runname
        from 
            run";
        
    for row in client.query(query, &[]).expect("Query failed.") {
            let info = RunInfo {
                runid: row.get(0), 
                runname: row.get(1)
            };
            println!("Run Info: {:?}", &info);
    }

}
