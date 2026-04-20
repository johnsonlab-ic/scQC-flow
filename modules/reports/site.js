document.addEventListener("DOMContentLoaded", () => {
  const body = document.body;
  const toggle = document.querySelector("[data-report-site-toggle]");
  const overlay = document.querySelector("[data-report-site-overlay]");
  const navLinks = document.querySelectorAll(".report-site-nav-link");

  function setNavOpen(isOpen) {
    body.classList.toggle("report-site-nav-open", isOpen);
    if (toggle) {
      toggle.setAttribute("aria-expanded", String(isOpen));
    }
    if (overlay) {
      overlay.hidden = !isOpen;
    }
  }

  if (toggle) {
    toggle.addEventListener("click", () => {
      setNavOpen(!body.classList.contains("report-site-nav-open"));
    });
  }

  if (overlay) {
    overlay.addEventListener("click", () => {
      setNavOpen(false);
    });
  }

  for (const link of navLinks) {
    link.addEventListener("click", () => {
      if (window.innerWidth <= 960) {
        setNavOpen(false);
      }
    });
  }

  window.addEventListener("resize", () => {
    if (window.innerWidth > 960) {
      setNavOpen(false);
    }
  });
});