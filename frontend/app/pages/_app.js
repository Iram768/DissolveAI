// frontend/frontend/pages/_app.js
import "../styles/globals.css";
import { useEffect } from "react";

function MyApp({ Component, pageProps }) {
  useEffect(() => {
    // Disable right click
    document.addEventListener("contextmenu", (e) => e.preventDefault());

    // Disable Ctrl+P (print)
    document.addEventListener("keydown", (e) => {
      if (e.ctrlKey && e.key === "p") {
        e.preventDefault();
        alert("Printing is disabled on this site.");
      }
    });
  }, []);

  return <Component {...pageProps} />;
}

export default MyApp;
