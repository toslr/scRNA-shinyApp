function scrollToSection(sectionId) {
  const element = document.getElementById(sectionId);
  const topBar = document.getElementById('topbar');
  const topBarHeight = topBar.offsetHeight;
  const elementPosition = element.getBoundingClientRect().top + window.pageYOffset;
  
  // Add 20px extra padding
  const offsetPosition = elementPosition - topBarHeight - 20;

  window.scrollTo({
    top: offsetPosition,
    behavior: 'smooth'
  });
}

$(document).on('change', '.condition-management-checkbox', function() {
  // Get the checkbox state
  var isChecked = $(this).prop('checked');
  
  // Get the input ID from the DOM
  var inputId = $(this).attr('id');
  
  // Delay the Shiny value update slightly to avoid race conditions
  setTimeout(function() {
    Shiny.setInputValue(inputId, isChecked);
  }, 10);
});

// Initialize collapsible sidebar sections when the document is ready
$(document).ready(function() {
  // Set up click handlers for collapsible sections
  $('.collapsible-section .section-header').on('click', function() {
    var section = $(this).parent();
    section.toggleClass('collapsed');
    
    // Store the state in localStorage to remember it between sessions
    var sectionId = section.index();
    localStorage.setItem('section_' + sectionId + '_collapsed', section.hasClass('collapsed'));
  });
  
  // Initialize sections based on saved state (if any)
  $('.collapsible-section').each(function(index) {
    var defaultCollapsed = [2,3,4,5,6];
    var isCollapsed = defaultCollapsed.includes(index);
    if (isCollapsed) {
      $(this).addClass('collapsed');
    }
  });
  
  // Existing scroll function from your original script.js
  window.scrollToSection = function(sectionId) {
    var section = document.getElementById(sectionId);
    if (section) {
      section.scrollIntoView({ behavior: 'smooth' });
    }
  };
  
  // Monitor for new condition checkboxes being added to the DOM
  const observer = new MutationObserver(function(mutations) {
    mutations.forEach(function(mutation) {
      if (mutation.addedNodes && mutation.addedNodes.length > 0) {
        // Check if any condition checkboxes were added
        const checkboxes = $(mutation.addedNodes).find('input[id^="conditionManagement-active_"]');
        if (checkboxes.length > 0) {
          // Add our special class to these checkboxes for the event handler
          checkboxes.addClass('condition-management-checkbox');
        }
      }
    });
  });
  
  // Start observing the document body for changes
  observer.observe(document.body, { 
    childList: true, 
    subtree: true 
  });
});

// Handler for updating sample checkboxes programmatically
Shiny.addCustomMessageHandler("updateSampleCheckboxes", function(message) {
  const selectedSamples = message.samples;
  
  // Find all sample checkboxes
  const checkboxes = document.querySelectorAll('input.sample-select');
  
  // Update each checkbox
  checkboxes.forEach(function(checkbox) {
    const gsm = checkbox.getAttribute('data-gsm');
    checkbox.checked = selectedSamples.includes(gsm);
  });
  
  // Also update the "select all" checkbox
  const selectAllCheckbox = document.getElementById('select-all-samples');
  if (selectAllCheckbox) {
    selectAllCheckbox.checked = checkboxes.length > 0 && 
      Array.from(checkboxes).every(cb => cb.checked);
  }
});