
use rpassword::prompt_password;
use postgres::{Client, NoTls};


fn main() {

    // For security reasons, don't hardcode the user ID and password

    let username = prompt_password("Enter server username: ").expect("Failed to read username.");
    let password = prompt_password("Enter server password: ").expect("Failed to read password.");

    // Connect to the server

    let connection_string = format!("postgresql://{}:{}@gaiadb08i:55431/surveys", username, password);
    let mut client = Client::connect(&connection_string, NoTls).expect("Failed to connect to the server.");


    println!("Hello, world!");

}
