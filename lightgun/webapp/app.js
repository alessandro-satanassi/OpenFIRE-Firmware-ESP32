document.addEventListener("DOMContentLoaded", () => {
    const statusText = document.getElementById("status");
    const btnTest = document.getElementById("btn-test");

    // Simulazione di connessione WebSocket
    statusText.innerText = "Loaded! Waiting for WebSocket...";
    statusText.style.color = "#4CAF50";

    btnTest.addEventListener("click", () => {
        alert("Pew pew! Test button clicked.");
    });
});
